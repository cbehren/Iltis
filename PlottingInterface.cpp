#include "PlottingInterface.H"
#include "SphericalShellData.H"
#include "InfiniteSlabData.H"
#include "RandomNumbers.H"
#include "LightParmParse/LParmParse.H"
#include "LymanAlphaLine.H"
#include "Parallel.H"
#include "Utilities.H"
#include "ListEmissionModel.H"
#include "UnigridDataset.H"
#include "Selectors.H"


int PlottingInterface::generate_rays()
{
    line_of_sight up2;
    normalize(los.x);
    normalize(up.x);
    crossproduct(los.x,up.x,up2.x);
    normalize(up2.x);
    //get the lower left corner of the field.
    double lowerleft[3];
    int lindim = npixels*oversampling;
    for(int i=0;i<3;i++)
    {
        lowerleft[i] = center.x[i] - up.x[i]*(width/2.0) - up2.x[i]*(width/2.0) - los.x[i]*(depth/2.0);
    }
    double dx = width/(npixels*oversampling);
    for(int i=0;i<npixels*oversampling;i++)
    {
        for(int j=0;j<npixels*oversampling;j++)
        {
            BaseParticle p;
            for(int dim=0;dim<3;dim++)
            {
                p.x[dim] =  lowerleft[dim];
                p.x[dim] += up.x[dim]*dx*i;
                p.x[dim] += up2.x[dim]*dx*j;
                p.lsp[dim] += up2.x[dim]*dx*j;
                p.k[dim] = los.x[dim];
                p.status = BP_OKAY;
                p.id = i*lindim+j;
                
            }
            
            
            //only IO proc stores the photons to avoid doubling.
            if(Parallel::IOProcessor())
            {
                //std::cout << p << std::endl;
                particles.push_back(p);
            }
            
        
        }
    
        
        
    }
    ds->set_domain(&particles);
    particles.exchange();
    
    
    return 0;
    
    
}



int PlottingInterface::do_step()
{
    do_rays();
    
    do_communication();
    if(verbosity>0)
    {
        long size = particles.size();
        Parallel::ReduceLongSum(size);
        if(Parallel::IOProcessor())
            std::cout << size << " photons left." << std::endl;
    }
    return 0;
}

bool PlottingInterface::is_done()
{
    
    bool value=false;
    if(particles.size()==0)
        value=true;
    //check if all processes are done.
    Parallel::ReduceBoolAnd(value);
            
    return value;
    
    
}
#define TRACK_STUCK_PHOTONS
#define RAY_DEBUG 0
//#define NO_SCATTERINGS
int PlottingInterface::do_rays()
{
    std::cout.precision(10);
    if(verbosity>1 && Parallel::IOProcessor())
        std::cout << "Entering do_rays" << std::endl;
    //TODO sort by order to reduce caching effort
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for(unsigned long i=0;i<particles.size();i++)
    {
        //main loop.
        BaseParticle &p = particles[i];

        DatasetStatus status;
        double pathlength=0.0;
        const BaseCell* cell = ds->data_at(p.x,p.k,pathlength,status);
        double dx = ds->get_dx(p.x,status);
#ifdef TRACK_STUCK_PHOTONS
        int num_steps_this_cell=0;  
        const BaseCell* old_cell = cell;
#endif
        bool go_on = true;
        while(status==DS_ALL_GOOD && go_on)
        {
            #ifdef TRACK_STUCK_PHOTONS
            if(old_cell==cell)
                num_steps_this_cell++;
            else
                num_steps_this_cell=0;
            if(num_steps_this_cell>1e5)
                std::cout << "many steps here" << p << " pathlength " << pathlength << std::endl;
            old_cell = cell;
            #endif
            //cycle throught the operators, store the result accordingly
            for(unsigned int iop=0;iop<plotops.size();iop++)
            {
                plotops[iop]->increment(p.id,&p,cell,pathlength,dx);
                //std::cout << dx << " " << output[iop][p.id] << " " << weight << std::endl;
            }
   
            if(RAY_DEBUG) 
                std::cout << "rho " << cell->density << " vel " << cell->velocity[0] << " " << cell->velocity[1]<< " " << cell->velocity[2] << " len " << pathlength <<  std::endl;
                
            p.move(pathlength);
            if(p.path_length>depth)
            {
                p.status = BP_DONE;
                go_on = false;
                break;
             }
            if(RAY_DEBUG) 
                    std::cout << "SCATTERED " << p.frequency << " #" << p.number_of_scatterings << " at " << p.x[0] << " " << p.x[1] << " " << p.x[2] << " k " << p.k[0] << " " << p.k[1] << " " << p.k[2] << " optical depth " << p.optical_depth_seen << std::endl;
                
            
            
            cell = ds->data_at(p.x,p.k,pathlength,status);
        }
        
        if(status == DS_NOT_IN_DOMAIN)
        {
            //TODO add periodic boundaries
            p.status = BP_DONE;
            
        }
        
    }
    return 0;
    
}

int PlottingInterface::do_communication()
{
    double start = Parallel::second();
    double v[3];
    if(verbosity>1 && Parallel::IOProcessor())
        std::cout << "Entering do_communication" << std::endl;
    //sort particles by their status. CAUTION: It is important to note that the order in which the ParticleStatus enum is defined is actually crucial here, because it moves the DEAD particles to the bottom, the DONE particles above them, and the OKAY particles to the top. 
    particles.sort_by_status();
    //do output for photons
    for(auto it = particles.rbegin(); it != particles.rend();it++)
    { 
     
        BaseParticle &p = *it;
        if(p.status == BP_DEAD)
        {
            //add hubble flow to frequency
            BaseEmissionLine::getLocalHubbleFlow(p,v);
            p.frequency = p.frequency - scalar(v,p.k)/Constants::CGS::cval*line->Nu0();
            particles.pop_back();
        }
        else if(p.status == BP_DONE)
        {
            //add hubble flow to frequency
            BaseEmissionLine::getLocalHubbleFlow(p,v);
            p.frequency = p.frequency - scalar(v,p.k)/Constants::CGS::cval*line->Nu0();
            particles.pop_back();
        }
    
    }

    //figure out which process should get which photon, and exchange photons accordingly
    ds->set_domain(&particles);
    {
        particles.exchange();
        if(verbosity>3)
            particles.print_balance();
        
    }
    if(verbosity > 1)
    {
        double end = Parallel::second();
        if(Parallel::IOProcessor())
        {
            std::cout << "Spent " << end-start << " seconds in communication " << std::endl;

        }

    }
    return 0;
}

int PlottingInterface::read_params()
{
    LParmParse pp;
    dataset_type = "SphericalShell";
    pp.query("dataset_type",dataset_type);
    pp.query("verbosity",verbosity);
    
    pp.get_line_of_sight("plotter.los",los);
    pp.get_line_of_sight("plotter.center",center);
    pp.get_line_of_sight("plotter.up",up);
    pp.get("plotter.width",width);
    pp.get("plotter.depth",depth);
    pp.get("plotter.npixels",npixels);
    pp.get("plotter.oversampling",oversampling);
    

    return 0;
}
int PlottingInterface::setup()
{
    LParmParse pp;

    ds = get_dataset(dataset_type);
    ds->setup();
    
    
    line = new LymanAlphaLine;
    line->setup();
    e = new BaseEmissionModel();
    o = new BaseOutput;
    if(Parallel::IOProcessor())
        std::cout << "PlottingInterface::setup(): Init complete. " << std::endl;
    
    pp.get("plotter.output_prefix",output_prefix);
    //generate the rays:
    generate_rays();
    //figure out what quantity we want to track.
    int k=0;
    std::string operator_name;
    while(pp.querykth("plotter.operator",k,operator_name))
    {
        std::cout << "Asked for operator " << operator_name << std::endl;
        
        plotops.push_back(get_plot_operator(operator_name,line));
        k++;
    }
    if(plotops.empty())
    {
        std::cout << "Error: No operator specified!" << std::endl;
        std::abort();
        
    }
    
    return 0;
}
int PlottingInterface::cleanup()
{
    
    for(unsigned int iop=0;iop<plotops.size();iop++)
    {
                plotops[iop]->write_results();
                //std::cout << dx << " " << output[iop][p.id] << " " << weight << std::endl;
    }
    return 0;
}

