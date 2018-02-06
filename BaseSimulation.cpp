#include "BaseSimulation.H"
#include "RandomNumbers.H"
#include "LightParmParse/LParmParse.H"
#include "LymanAlphaLine.H"
#include "Parallel.H"
#include "Utilities.H"
#include "ListEmissionModel.H"
#include "UnigridDataset.H"
#include "Selectors.H"





int BaseSimulation::do_step()
{
    do_rays();
    do_peeling_off();
    
    do_communication();
    if(verbosity>0)
    {
        long size = particles.size();
        Parallel::ReduceLongSum(size);
        if(Parallel::IOProcessor())
            std::cout << size << " photons left." << std::endl;
        if(use_peeling_off)
        {
            size = 0;
            for(int i=0;i<number_of_instruments;i++)
                size += vparticles[i].size();
             Parallel::ReduceLongSum(size);
            if(Parallel::IOProcessor())
                std::cout << size << " peeling off photons left." << std::endl;
        }
    }
    return 0;
}

bool BaseSimulation::is_done()
{
    
    bool value=false;
    if(particles.size()==0)
        value=true;
    if(value && use_peeling_off)
    {
        for(int i=0;i<number_of_instruments;i++)
            if(vparticles[i].size()==0)
                value=true;
            else
            {
                value=false;
                break;
            }
    }
    //check if all processes are done.
    Parallel::ReduceBoolAnd(value);
            
    return value;
    
    
}
#define TRACK_STUCK_PHOTONS
#define RAY_DEBUG 0
//#define NO_SCATTERINGS
int BaseSimulation::do_rays()
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
        double fraction_in_dust=0.0;
        const BaseCell* cell = ds->data_at(p.x,p.k,pathlength,status);
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
            if(p.number_of_scatterings == -1)
            {
                //this photon is new (signaled by the -1) - launch peeling off from its emission spot.
                launch_peeling(p,*cell,NULL);
                p.number_of_scatterings = 0;
                p.optical_depth = RNG::exponential();
                if(RAY_DEBUG) 
                    std::cout  <<  "EMITTED " << p.frequency << " #" << p.number_of_scatterings << " at " << p.x[0] << " " << p.x[1] << " " << p.x[2] << " k " << p.k[0] << " " << p.k[1] << " " << p.k[2] << " optical depth " << p.optical_depth << std::endl;

                
            }
            //get the optical depth, figure out whether we will make it through pathlength.
            double dtau0= line->get_optical_depth_line_center(cell,&p,pathlength,&fraction_in_dust);
            double dtau = line->get_optical_depth(cell,&p,pathlength,&fraction_in_dust);
//             if(RAY_DEBUG) std::cout << "photon " << p <<" in cell " << cell << " " << *cell << " mytau " << p.optical_depth << std::endl;
#ifdef REMOVE_CGM            
            //TODO just a try, remove
            double r=0.;
            for(int j=0;j<3;j++)
            {
                r+=pow(center[j]-p.x[j],2);
                
            }
            if(sqrt(r) > radius)
            {
                dtau0 = 0;
            }
            //TODO
#endif            
            if(RAY_DEBUG) 
            std::cout << "rho " << cell->density << " vel " << cell->velocity[0] << " " << cell->velocity[1]<< " " << cell->velocity[2] << " len " << pathlength << " dtau " << dtau <<  std::endl;
            
#ifndef NO_SCATTERINGS
            if(p.optical_depth >= dtau)
#else
            if(1)
#endif
            {
                
                //we make it - move the photon to its new position, decrease the remaining optical depth accordingly.
                p.move(pathlength);
                p.optical_depth -= dtau;
                p.optical_depth_seen += dtau;
//                 if(RAY_DEBUG) std::cout << " goes " << pathlength << " total " << p.path_length << "\tdtau " << dtau << "\ttauseen " << p.optical_depth_seen <<  std::endl;
            }
            else
            {
                //we will not make it, an interaction will occur.
#ifdef TRACK_STUCK_PHOTONS
                num_steps_this_cell=0;  
        
#endif                
                //move us to the point of interaction
                double pathlength2 = p.optical_depth/dtau*pathlength;
                if(RAY_DEBUG) 
                std::cout << " real len " << pathlength2<< " taur " << p.optical_depth << std::endl;
                p.move(pathlength2);
                p.optical_depth_seen += p.optical_depth;
                //execute the scattering
                double new_k[3];
                double atom_velocity[3];
#ifndef NO_SCATTERINGS
                double new_freq_obs;
                //decide whether we scatter on dust or gas
                double r = RNG::uniform();
                if(r <= fraction_in_dust)//we interact with dust!
                {
                    double r2 = RNG::uniform();
                    if(r2 > line->dm->getAlbedo())
                    {
                        //we are absorbed
                        if(use_biasing && p.bias>minimum_bias)
                        {
                            p.bias*=line->dm->getAlbedo();
                            new_freq_obs = line->dm->scatter(cell,&p,new_k,false);     
                            launch_peeling(p,*cell,NULL,true);
                        }
                        else
                        {
                            p.status = BP_DEAD;
                            go_on = false;
                            break;
                        }
                    }
                    else
                    {
                        //we are scattered on dust
                        new_freq_obs = line->dm->scatter(cell,&p,new_k,false);     
                        launch_peeling(p,*cell,NULL,true);
                    }
                    
                    
                }
                else //we interact with hydrogen!
                {                               
                    new_freq_obs = line->scatter(cell,&p,new_k,dtau0,atom_velocity,false);
                    launch_peeling(p,*cell,atom_velocity);
                }
                //update particle
                for(int j=0;j<3;j++)
                {
                    p.k[j] = new_k[j];
                    p.lsp[j] = p.x[j];
                    
                }
                p.number_of_scatterings++;
                p.frequency = new_freq_obs;
#endif //NO_SCATTERINGS
                p.optical_depth = RNG::exponential();
                if(RAY_DEBUG) 
                    std::cout << "SCATTERED " << p.frequency << " #" << p.number_of_scatterings << " at " << p.x[0] << " " << p.x[1] << " " << p.x[2] << " k " << p.k[0] << " " << p.k[1] << " " << p.k[2] << " optical depth " << p.optical_depth_seen << std::endl;
                
            }
            if(use_peeling_off && vparticles[0].size()*number_of_instruments > (unsigned long)max_num_peeling_off_photons)
                break;
            
            cell = ds->data_at(p.x,p.k,pathlength,status);
            
        }
        if(status == DS_NOT_IN_DOMAIN)
        {
            //TODO add periodic boundaries
            p.status = BP_DONE;
            
        }
        //if we have too many peeling off photons, stop propagating.
        if(use_peeling_off && vparticles[0].size() > (unsigned long)max_num_peeling_off_photons)
                continue;
        
    }
    return 0;
    
}
int BaseSimulation::do_peeling_off()
{
    if(verbosity>1 && Parallel::IOProcessor())
        std::cout << "Entering do_peeling_off" << std::endl;
        //TODO sort by order to reduce caching effort
    double start = Parallel::second();
    long size = 0;
    //cycle over instruments
    for(int j=0;j<number_of_instruments;j++)
    {
        BaseParticleVector& vparticle = vparticles[j];
        if(verbosity>2)
        {
            size = vparticle.size();
            Parallel::ReduceLongSum(size);
            if(Parallel::IOProcessor())
                std::cout << size << " peeling off photons in queue for instrument " << j << std::endl;
            
        }
        //cycle over peeling off photons.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for(unsigned long i=0;i<vparticle.size();i++)
        {
            //main loop.
            BaseParticle &p = vparticle[i];
            
            DatasetStatus status;
            double pathlength=0.0;
            double fraction_in_dust=0.0;
            const BaseCell* cell = ds->data_at(p.x,p.k,pathlength,status);
            bool go_on = true;
            while(status==DS_ALL_GOOD && go_on)
            {
                double dtau = line->get_optical_depth(cell,&p,pathlength,&fraction_in_dust);
#ifdef REMOVE_CGM            
            //TODO just a try, remove
            double r=0.;
            for(int j=0;j<3;j++)
            {
                r+=pow(center[j]-p.x[j],2);
                
            }
            if(sqrt(r) > radius)
            {
                dtau = 0;
            }
#endif            
                p.move(pathlength);
                //accumulate optical depth seen.
                p.optical_depth += dtau;
                p.optical_depth_seen += dtau;
                cell = ds->data_at(p.x,p.k,pathlength,status);
                //discard photon is tau becomes too large.
                if(p.optical_depth > tau_max)
                {
                    p.status = BP_DEAD;
                    go_on = false;
                }

            }
            if(status == DS_NOT_IN_DOMAIN)
            {
                //TODO periodic boundaries
                p.status = BP_DONE;
                
            }
        }
    }
    Parallel::Barrier();
    if(verbosity > 1)
    {
        double end = Parallel::second();
        if(Parallel::IOProcessor())
        {
            std::cout << "Spent " << end-start << " seconds in peeling off, " << size*number_of_instruments/(end-start) << " photons per second." << std::endl;

        }

    }
    return 0;
    
}

int BaseSimulation::do_communication()
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
            o->writeDeadPhoton(p);
            particles.pop_back();
        }
        else if(p.status == BP_DONE)
        {
            //add hubble flow to frequency
            BaseEmissionLine::getLocalHubbleFlow(p,v);
            p.frequency = p.frequency - scalar(v,p.k)/Constants::CGS::cval*line->Nu0();
            o->writePhoton(p);
            particles.pop_back();
        }
    
    }
    //do output for peeling off photons
    for(int j=0;j<number_of_instruments;j++)
    {
        vparticles[j].sort_by_status();
        for(auto it = vparticles[j].rbegin(); it != vparticles[j].rend();it++)
        { 
        
            BaseParticle &p = *it;
            if(p.status == BP_DEAD)
            {
                vparticles[j].pop_back();
            }
            else if(p.status == BP_DONE)
            {
                //add hubble flow to frequency
                BaseEmissionLine::getLocalHubbleFlow(p,v);
                p.frequency = p.frequency - scalar(v,p.k)/Constants::CGS::cval*line->Nu0();
                o->writePeelingOff(p,j);
                vparticles[j].pop_back();
            }
        
        }
    }
    //figure out which process should get which photon, and exchange photons accordingly
    ds->set_domain(&particles);
    {
        particles.exchange();
        if(verbosity>3)
            particles.print_balance();
        
    }
    //same for all the instruments for the peeling off photons.
    for(int j=0;j<number_of_instruments;j++)
    {
        ds->set_domain(&vparticles[j]);
        vparticles[j].exchange(); 
        if(Parallel::IOProcessor() && verbosity>3)
            std::cout << "For instrument " << j << ": " << std::endl;
        vparticles[j].print_balance();
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

int BaseSimulation::read_params()
{
    LParmParse pp;
    dataset_type = "SphericalShell";
    pp.query("dataset_type",dataset_type);
    pp.query("use_peeling_off",use_peeling_off);
    pp.query("verbosity",verbosity);
    pp.query("max_num_peeling_off_photons",max_num_peeling_off_photons);
    std::string dummy;
    if(pp.query("emission.emitter_file",dummy))
        use_emitter_file = true;
    return 0;
}
int BaseSimulation::setup()
{
    LParmParse pp;
    if(use_peeling_off)
    {
        number_of_instruments = 1;
        pp.query("number_of_instruments",number_of_instruments);
        vparticles = new BaseParticleVector[number_of_instruments];
        std::vector<double> los;

        for(int i=0;i<number_of_instruments;i++)
        {
            pp.getktharr("line_of_sight",i,los,LParmParse::FIRST,3);
            if(los.size()!=3)
            {
                std::cout << "Failed to read instrument." << std::endl;
                std::abort();
                
            }
            line_of_sight mylos;
            for(int j=0;j<3;j++)
            {
                mylos.x[j] = los[j];
            }
            normalize(mylos.x);
            if(Parallel::IOProcessor())
                std::cout << "read instrument " << i << " with line of sight (" << mylos.x[0] << " " << mylos.x[1] << " " << mylos.x[2] << ") " << std::endl;
            instruments.push_back(mylos);
        }
        if(pp.queryktharr("line_of_sight",number_of_instruments,los))
        {
            std::cout << "Too many LoS in inputs file." << std::endl;
            std::abort();
        }
        tau_max = 20.0;
        pp.query("tau_max",tau_max);
        
    }
    
    pp.query("use_biasing",use_biasing);
    if(use_biasing)
        if(pp.query("minimum_bias",minimum_bias))
            if(Parallel::IOProcessor())
                std::cout << "Using minimum bias of " << minimum_bias << std::endl;
    
    ds = get_dataset(dataset_type);
    ds->setup();
    
    line = new LymanAlphaLine;
    line->setup();
    
    
    if(use_emitter_file)
    {
        e = new ListEmissionModel();
        if(Parallel::IOProcessor())
            std::cout << "Using emitters from file" << std::endl;
    }
    else
        e = new BaseEmissionModel();
    e->setup(particles,*line,*ds);
    
    o = new BaseOutput;
    o->setup(line);
    o->writeSlices(*ds);
    o->writeInputPhotons(particles);
    
    if(Parallel::IOProcessor())
        std::cout << "BaseSimulation::setup(): Init complete. " << std::endl;
    
    return 0;
}
int BaseSimulation::cleanup()
{
    o->cleanup();
    return 0;
}
void BaseSimulation::launch_peeling(const BaseParticle &p, const BaseCell& cell, double atom_velocity[3],bool on_dust)
{
    for(int i=0;i<number_of_instruments;i++)
    {
        BaseParticle vp=p;
        for(int j=0;j<3;j++)
        {
            vp.k[j] = instruments[i].x[j];
            vp.lsp[j] = vp.x[j];
        }
        vp.optical_depth = 0.0;
        vp.optical_depth_seen = 0.0;
        
        double new_freq_obs=0.0;
        if(on_dust)
        {
            new_freq_obs = line->dm->scatter(&cell,&p,vp.k,true);
            //TODO to we need to multiply by the albedo or not?
            vp.weight *= line->dm->GreensteinProbability(scalar(vp.k,p.k));
        }
        else
        {
            if(atom_velocity)
                new_freq_obs = line->scatter(&cell,&p,vp.k,0.0,atom_velocity,true);
            else //the photon in question has just been emitted; we need to shift around in frequency.
            {
                new_freq_obs = p.frequency +  (-scalar(cell.velocity,p.k)/ Constants::CGS::cval + scalar(cell.velocity,vp.k)/Constants::CGS::cval)*line->Nu0();
                
            }
        }
        vp.frequency = new_freq_obs;
        double pathlength = 0.0;

        //we can save us some trouble here.
        double fraction_in_dust = 0.0;
        double dtau = line->get_optical_depth(&cell,&vp,pathlength,&fraction_in_dust);
        if(dtau>tau_max)
            continue;
        //std::cout << vp;
#pragma omp critical 
        {
            vparticles[i].push_back(vp);
        }
    }
}



