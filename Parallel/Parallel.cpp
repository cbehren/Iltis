
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <ctime>
#include "sys/time.h"

#include "Parallel.H"

//This is a stripped/modified version of the ParallelDescriptor module as included in the BoxLib library, see https://github.com/BoxLib-Codes/BoxLib 
//For details of the license, see license.txt

namespace Parallel
{
    //
    // My processor ID.
    //
    int m_MyId = -1;
    //
    // The number of processors.
    //
    int m_nProcs = -1;
    //
    // Communicator
    //
    MPI_Comm m_comm;

    const int ioProcessor = 0;

    namespace util
    {
	//
	// Reduce helper functons.
	//
	void DoAllReduceDouble (double& r, MPI_Op op);
	void DoAllReduceLong (long& r, MPI_Op op);
	void DoAllReduceInt  (int&  r, MPI_Op op);

	void DoAllReduceDouble (double* r, MPI_Op op, int cnt);
	void DoAllReduceLong (long* r, MPI_Op op, int cnt);
	void DoAllReduceInt  (int*  r, MPI_Op op, int cnt);

	void DoReduceDouble (double& r, MPI_Op op, int cpu);
	void DoReduceLong (long& r, MPI_Op op, int cpu);
	void DoReduceInt  (int&  r, MPI_Op op, int cpu);

	void DoReduceDouble (double* r, MPI_Op op, int cnt, int cpu);
	void DoReduceLong (long* r, MPI_Op op, int cnt, int cpu);
	void DoReduceInt  (int*  r, MPI_Op op, int cnt, int cpu);
    }
}

//
// Definition of non-inline members of CommData.
//


#ifdef USE_MPI

#include <mpi.h>

namespace
{
    const char*
    the_message_string (const char* file,
                        int         line,
                        const char* call,
                        int         status)
    {
	const int N = 512;
	static char buf[N];
	if ( status )
	{
	    snprintf(buf, N, "MPI Error: File %s, line %d, %s: %s",
                     file, line, call, Parallel::ErrorString(status));
	}
	else
	{
	    snprintf(buf, N, "MPI Error: File %s, line %d, %s",
                     file, line, call);
	}
	return buf;
    }

}

namespace Parallel
{
    void
    MPI_Error (const char* file, int line, const char* str, int rc)
    {
        std::cout << (the_message_string(file, line, str, rc));
        std::abort();
    }
}

void
Parallel::Abort ()
{
    MPI_Abort(Communicator(), -1);
}

void
Parallel::Abort (int errorcode)
{

    MPI_Abort(Communicator(), errorcode);

}

const char*
Parallel::ErrorString (int errorcode)
{

    int len = 0;

    static char msg[MPI_MAX_ERROR_STRING+1];

    MPI_Error_string(errorcode, msg, &len);

    return msg;
}

void
Parallel::Message::wait ()
{

    BL_MPI_REQUIRE( MPI_Wait(&m_req, &m_stat) );
}

bool
Parallel::Message::test ()
{
    int flag;
    BL_MPI_REQUIRE( MPI_Test(&m_req, &flag, &m_stat) );
    m_finished = flag != 0;
    return m_finished;
}

int
Parallel::Message::tag () const
{
    if ( !m_finished ) 
    {
     std::cout <<  "Message::tag: Not Finished!";
     std::abort();
    }
    return m_stat.MPI_TAG;
}

int
Parallel::Message::pid () const
{
    if ( !m_finished ) 
    {
        std::cout << "Message::pid: Not Finished!";
        std::abort();
    }
    return m_stat.MPI_SOURCE;
}

size_t
Parallel::Message::count () const
{
    if ( m_type == MPI_DATATYPE_NULL ) 
    {
        std::cout << "Message::count: Bad Type!";
        std::abort();
    }
    if ( !m_finished ) 
    {
       std::cout << "Message::count: Not Finished!";
       std::abort();
    }
    int cnt;
    BL_MPI_REQUIRE( MPI_Get_count(&m_stat, m_type, &cnt) );
    return cnt;
}
//CHANGE by Christoph Behrens: Add additional optional argument to start with threads.
void
Parallel::StartParallel (int*    argc,
                                   char*** argv,
                                   MPI_Comm mpi_comm,
                                   int mode_required
                                  )
{
    m_comm = mpi_comm;
    int sflag;

    BL_MPI_REQUIRE( MPI_Initialized(&sflag) );
    int provided=-1;

    if (!sflag)
	BL_MPI_REQUIRE( MPI_Init_thread(argc, argv, mode_required, &provided) );
    if(mode_required != provided)
    {
        std::cout << "StartParallel: Could not set the required level of MPI_THREAD_* \n";
        std::abort();
    }
    BL_MPI_REQUIRE( MPI_Comm_size(Communicator(), &m_nProcs) );

    BL_MPI_REQUIRE( MPI_Comm_rank(Communicator(), &m_MyId) );
    //
    // Wait till all other processes are properly started.
    //
    BL_MPI_REQUIRE( MPI_Barrier(Communicator()) );
}

void
Parallel::EndParallel ()
{
    
    BL_MPI_REQUIRE( MPI_Finalize() );
}

double
Parallel::second ()
{
    return MPI_Wtime();
}

void
Parallel::Barrier ()
{
    
    BL_MPI_REQUIRE( MPI_Barrier(Parallel::Communicator()) );
}

void
Parallel::Barrier (MPI_Comm comm)
{
    
    BL_MPI_REQUIRE( MPI_Barrier(comm) );
}

void
Parallel::Test (MPI_Request& request, int& flag, MPI_Status& status)
{
    BL_MPI_REQUIRE( MPI_Test(&request,&flag,&status) );
}

void
Parallel::IProbe (int src_pid, int tag, int& flag, MPI_Status& status)
{
    BL_MPI_REQUIRE( MPI_Iprobe(src_pid, tag, Parallel::Communicator(),
                               &flag, &status) );
}

void
Parallel::Comm_dup (MPI_Comm comm, MPI_Comm& newcomm)
{
    BL_MPI_REQUIRE( MPI_Comm_dup(comm, &newcomm) );
}

void
Parallel::util::DoAllReduceDouble (double&  r,
                                           MPI_Op op)
{
   
    double recv;

    BL_MPI_REQUIRE( MPI_Allreduce(&r,
                                  &recv,
                                  1,
                                  Mpi_typemap<double>::type(),
                                  op,
                                  Communicator()) );
    r = recv;
}

void
Parallel::util::DoAllReduceDouble (double*  r,
                                           MPI_Op op,
                                           int    cnt)
{


    std::vector<double> recv(cnt);

    BL_MPI_REQUIRE( MPI_Allreduce(r,
                                  recv.data(),
                                  cnt,
                                  Mpi_typemap<double>::type(),
                                  op,
                                  Communicator()) );
    for (int i = 0; i < cnt; i++)
        r[i] = recv[i];
}

void
Parallel::util::DoReduceDouble (double&  r,
                                        MPI_Op op,
                                        int    cpu)
{

    double recv;

    BL_MPI_REQUIRE( MPI_Reduce(&r,
                               &recv,
                               1,
                               Mpi_typemap<double>::type(),
                               op,
                               cpu,
                               Communicator()) );

    if (Parallel::MyProc() == cpu)
        r = recv;
}

void
Parallel::util::DoReduceDouble (double*  r,
                                        MPI_Op op,
                                        int    cnt,
                                        int    cpu)
{
    std::vector<double> recv(cnt);

    BL_MPI_REQUIRE( MPI_Reduce(r,
                               recv.data(),
                               cnt,
                               Mpi_typemap<double>::type(),
                               op,
                               cpu,
                               Communicator()) );

    if (Parallel::MyProc() == cpu)
    {
        for (int i = 0; i < cnt; i++)
            r[i] = recv[i];
    }
}

void
Parallel::ReduceDoubleMax (double& r)
{
    util::DoAllReduceDouble(r,MPI_MAX);
}

void
Parallel::ReduceDoubleMin (double& r)
{
    util::DoAllReduceDouble(r,MPI_MIN);
}

void
Parallel::ReduceDoubleSum (double& r)
{
    util::DoAllReduceDouble(r,MPI_SUM);
}

void
Parallel::ReduceDoubleMax (double* r, int cnt)
{
    util::DoAllReduceDouble(r,MPI_MAX,cnt);
}

void
Parallel::ReduceDoubleMin (double* r, int cnt)
{
    util::DoAllReduceDouble(r,MPI_MIN,cnt);
}

void
Parallel::ReduceDoubleSum (double* r, int cnt)
{
    util::DoAllReduceDouble(r,MPI_SUM,cnt);
}

void
Parallel::ReduceDoubleMax (double& r, int cpu)
{
    util::DoReduceDouble(r,MPI_MAX,cpu);
}

void
Parallel::ReduceDoubleMin (double& r, int cpu)
{
    util::DoReduceDouble(r,MPI_MIN,cpu);
}

void
Parallel::ReduceDoubleSum (double& r, int cpu)
{
    util::DoReduceDouble(r,MPI_SUM,cpu);
}

void
Parallel::ReduceDoubleMax (double* r, int cnt, int cpu)
{
    util::DoReduceDouble(r,MPI_MAX,cnt,cpu);
}


void
Parallel::ReduceDoubleMin (double* r, int cnt, int cpu)
{
    util::DoReduceDouble(r,MPI_MIN,cnt,cpu);
}

void
Parallel::ReduceDoubleSum (double* r, int cnt, int cpu)
{
    util::DoReduceDouble(r,MPI_SUM,cnt,cpu);
}

void
Parallel::util::DoAllReduceLong (long&  r,
                                           MPI_Op op)
{
    
    long recv;

    BL_MPI_REQUIRE( MPI_Allreduce(&r,
                                  &recv,
                                  1,
                                  MPI_LONG,
                                  op,
                                  Communicator()) );
    r = recv;
}

void
Parallel::util::DoAllReduceLong (long*  r,
                                           MPI_Op op,
                                           int    cnt)
{

    std::vector<long> recv(cnt);

    BL_MPI_REQUIRE( MPI_Allreduce(r,
                                  recv.data(),
                                  cnt,
                                  MPI_LONG,
                                  op,
                                  Communicator()) );
    for (int i = 0; i < cnt; i++)
        r[i] = recv[i];
}

void
Parallel::util::DoReduceLong (long&  r,
                                        MPI_Op op,
                                        int    cpu)
{
    long recv;

    BL_MPI_REQUIRE( MPI_Reduce(&r,
                               &recv,
                               1,
                               MPI_LONG,
                               op,
                               cpu,
                               Communicator()));

    if (Parallel::MyProc() == cpu)
        r = recv;
}

void
Parallel::util::DoReduceLong (long*  r,
                                        MPI_Op op,
                                        int    cnt,
                                        int    cpu)
{
    std::vector<long> recv(cnt);

    BL_MPI_REQUIRE( MPI_Reduce(r,
                               recv.data(),
                               cnt,
                               MPI_LONG,
                               op,
                               cpu,
                               Communicator()));

    if (Parallel::MyProc() == cpu)
    {
        for (int i = 0; i < cnt; i++)
            r[i] = recv[i];
    }
}

void
Parallel::ReduceLongAnd (long& r)
{
    util::DoAllReduceLong(r,MPI_LAND);
}

void
Parallel::ReduceLongSum (long& r)
{
    util::DoAllReduceLong(r,MPI_SUM);
}

void
Parallel::ReduceLongMax (long& r)
{
    util::DoAllReduceLong(r,MPI_MAX);
}

void
Parallel::ReduceLongMin (long& r)
{
    util::DoAllReduceLong(r,MPI_MIN);
}

void
Parallel::ReduceLongAnd (long* r, int cnt)
{
    util::DoAllReduceLong(r,MPI_LAND,cnt);
}

void
Parallel::ReduceLongSum (long* r, int cnt)
{
    util::DoAllReduceLong(r,MPI_SUM,cnt);
}

void
Parallel::ReduceLongMax (long* r, int cnt)
{
    util::DoAllReduceLong(r,MPI_MAX,cnt);
}

void
Parallel::ReduceLongMin (long* r, int cnt)
{
    util::DoAllReduceLong(r,MPI_MIN,cnt);
}

void
Parallel::ReduceLongAnd (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_LAND,cpu);
}

void
Parallel::ReduceLongSum (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_SUM,cpu);
}

void
Parallel::ReduceLongMax (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_MAX,cpu);
}

void
Parallel::ReduceLongMin (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_MIN,cpu);
}

void
Parallel::ReduceLongAnd (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_LAND,cnt,cpu);
}

void
Parallel::ReduceLongSum (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_SUM,cnt,cpu);
}

void
Parallel::ReduceLongMax (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_MAX,cnt,cpu);
}

void
Parallel::ReduceLongMin (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_MIN,cnt,cpu);
}

void
Parallel::util::DoAllReduceInt (int&   r,
                                          MPI_Op op)
{
    int recv;

    BL_MPI_REQUIRE( MPI_Allreduce(&r,
                                  &recv,
                                  1,
                                  MPI_INT,
                                  op,
                                  Communicator()));
    r = recv;
}

void
Parallel::util::DoAllReduceInt (int*   r,
                                          MPI_Op op,
                                          int    cnt)
{
    std::vector<int> recv(cnt);

    BL_MPI_REQUIRE( MPI_Allreduce(r,
                                  recv.data(),
                                  cnt,
                                  MPI_INT,
                                  op,
                                  Communicator()));
    for (int i = 0; i < cnt; i++)
        r[i] = recv[i];
}

void
Parallel::util::DoReduceInt (int&   r,
                                       MPI_Op op,
                                       int    cpu)
{
    int recv;

    BL_MPI_REQUIRE( MPI_Reduce(&r,
                               &recv,
                               1,
                               MPI_INT,
                               op,
                               cpu,
                               Communicator()));

    if (Parallel::MyProc() == cpu)
        r = recv;
}

void
Parallel::util::DoReduceInt (int*   r,
                                       MPI_Op op,
                                       int    cnt,
                                       int    cpu)
{

    std::vector<int> recv(cnt);

    BL_MPI_REQUIRE( MPI_Reduce(r,
                               recv.data(),
                               cnt,
                               MPI_INT,
                               op,
                               cpu,
                               Communicator()));

    if (Parallel::MyProc() == cpu)
    {
        for (int i = 0; i < cnt; i++)
            r[i] = recv[i];
    }
}

void
Parallel::ReduceIntSum (int& r)
{
    util::DoAllReduceInt(r,MPI_SUM);
}

void
Parallel::ReduceIntMax (int& r)
{
    util::DoAllReduceInt(r,MPI_MAX);
}

void
Parallel::ReduceIntMin (int& r)
{
    util::DoAllReduceInt(r,MPI_MIN);
}

void
Parallel::ReduceIntSum (int* r, int cnt)
{
    util::DoAllReduceInt(r,MPI_SUM,cnt);
}

void
Parallel::ReduceIntMax (int* r, int cnt)
{
    util::DoAllReduceInt(r,MPI_MAX,cnt);
}

void
Parallel::ReduceIntMin (int* r, int cnt)
{
    util::DoAllReduceInt(r,MPI_MIN,cnt);
}

void
Parallel::ReduceIntSum (int& r, int cpu)
{
    util::DoReduceInt(r,MPI_SUM,cpu);
}

void
Parallel::ReduceIntMax (int& r, int cpu)
{
    util::DoReduceInt(r,MPI_MAX,cpu);
}

void
Parallel::ReduceIntMin (int& r, int cpu)
{
    util::DoReduceInt(r,MPI_MIN,cpu);
}

void
Parallel::ReduceIntSum (int* r, int cnt, int cpu)
{
    util::DoReduceInt(r,MPI_SUM,cnt,cpu);
}

void
Parallel::ReduceIntMax (int* r, int cnt, int cpu)
{
    util::DoReduceInt(r,MPI_MAX,cnt,cpu);
}

void
Parallel::ReduceIntMin (int* r, int cnt, int cpu)
{
    util::DoReduceInt(r,MPI_MIN,cnt,cpu);
}

void
Parallel::ReduceBoolAnd (bool& r)
{
    int src = r; // src is either 0 or 1.

    util::DoAllReduceInt(src,MPI_SUM);

    r = (src == Parallel::NProcs()) ? true : false;
}

void
Parallel::ReduceBoolOr (bool& r)
{
    int src = r; // src is either 0 or 1.

    util::DoAllReduceInt(src,MPI_SUM);

    r = (src == 0) ? false : true;
}

void
Parallel::ReduceBoolAnd (bool& r, int cpu)
{
    int src = r; // src is either 0 or 1.

    util::DoReduceInt(src,MPI_SUM,cpu);

    if (Parallel::MyProc() == cpu)
        r = (src == Parallel::NProcs()) ? true : false;
}

void
Parallel::ReduceBoolOr (bool& r, int cpu)
{
    int src = r; // src is either 0 or 1.

    util::DoReduceInt(src,MPI_SUM,cpu);

    if (Parallel::MyProc() == cpu)
        r = (src == 0) ? false : true;
}


void Parallel::Gatherv (double* sendbuf,
                  int sendcount,
                  int*   recvcounts,
                  double* recvbuf,
                  int   root)
{
    MPI_Datatype typ = Mpi_typemap<double>::type();
    std::vector<int> offsets(NProcs());
    for(int i=1;i<NProcs();i++)
        offsets[i] = offsets[i-1]+recvcounts[i-1];
    BL_MPI_REQUIRE( MPI_Gatherv(sendbuf,
                               sendcount,
                               typ,
                               recvbuf,
                               recvcounts,
                                offsets.data(),
                               typ,
                               root,
                               Communicator()));
    
    
}

void Parallel::Gatherv (int* sendbuf,
                  int sendcount,
                  int*   recvcounts,
                  int* recvbuf,
                  int   root)
{
    MPI_Datatype typ = Mpi_typemap<int>::type();
    std::vector<int> offsets(NProcs());
    for(int i=1;i<NProcs();i++)
        offsets[i] = offsets[i-1]+recvcounts[i-1];
    BL_MPI_REQUIRE( MPI_Gatherv(sendbuf,
                               sendcount,
                               typ,
                               recvbuf,
                               recvcounts,
                                offsets.data(),
                               typ,
                               root,
                               Communicator()));
    
    
}

void Parallel::Gatherv (long* sendbuf,
                  int sendcount,
                  int*   recvcounts,
                  long* recvbuf,
                  int   root)
{
    MPI_Datatype typ = Mpi_typemap<long>::type();
    std::vector<int> offsets(NProcs());
    for(int i=1;i<NProcs();i++)
        offsets[i] = offsets[i-1]+recvcounts[i-1];
    BL_MPI_REQUIRE( MPI_Gatherv(sendbuf,
                               sendcount,
                               typ,
                               recvbuf,
                               recvcounts,
                                offsets.data(),
                               typ,
                               root,
                               Communicator()));
    
    
}

void
Parallel::Gather (double* sendbuf,
                            int   nsend,
                            double* recvbuf,
                            int   root)
{
    MPI_Datatype typ = Mpi_typemap<double>::type();

    BL_MPI_REQUIRE( MPI_Gather(sendbuf,
                               nsend,
                               typ,
                               recvbuf,
                               nsend,
                               typ,
                               root,
                               Communicator()));
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<char>::type ()
{
    return  MPI_CHAR;
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<short>::type ()
{
    return  MPI_SHORT;
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<int>::type ()
{
    return  MPI_INT;
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<long>::type ()
{
    return  MPI_LONG;
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<unsigned char>::type ()
{
    return  MPI_UNSIGNED_CHAR;
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<unsigned short>::type ()
{
    return  MPI_UNSIGNED_SHORT;
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<unsigned int>::type ()
{
    return  MPI_UNSIGNED;
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<unsigned long>::type ()
{
    return  MPI_UNSIGNED_LONG;
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<float>::type ()
{
    return  MPI_FLOAT;
}

template <>
MPI_Datatype
Parallel::Mpi_typemap<double>::type ()
{
    return  MPI_DOUBLE;
}

#else /*!USE_MPI*/


void Parallel::Gatherv (double* sendbuf,
                  int sendcount,
                  int*   recvcounts,
                  double* recvbuf,
                  int   root)
{
    
    for (int i = 0; i < sendcount; ++i)
        recvbuf[i] = sendbuf[i];
}

void Parallel::Gatherv (int* sendbuf,
                  int sendcount,
                  int*   recvcounts,
                  int* recvbuf,
                  int   root)
{
    
    for (int i = 0; i < sendcount; ++i)
        recvbuf[i] = sendbuf[i];
}

void Parallel::Gatherv (long* sendbuf,
                  int sendcount,
                  int*   recvcounts,
                  long* recvbuf,
                  int   root)
{
    
    for (int i = 0; i < sendcount; ++i)
        recvbuf[i] = sendbuf[i];
}
void
Parallel::StartParallel (int*    argc,
                                   char*** argv,
                                   MPI_Comm, int required)
{
    m_nProcs    = 1;
    m_MyId      = 0;
    m_comm      = 0;
}

void
Parallel::Gather (double* sendbuf,
			    int   nsend,
			    double* recvbuf,
			    int   root)
{
    for (int i = 0; i < nsend; ++i)
        recvbuf[i] = sendbuf[i];
}

void
Parallel::Message::wait ()
{}

bool
Parallel::Message::test ()
{
    return m_finished;
}

void Parallel::EndParallel () {}

void Parallel::Abort ()
{
    std::abort(); 
}
void Parallel::Abort (int)
{ 
    std::abort(); 
}

const char* Parallel::ErrorString (int) { return ""; }

void Parallel::Barrier () {}
void Parallel::Barrier (MPI_Comm) {}

void Parallel::Test (MPI_Request&, int&, MPI_Status&) {}
void Parallel::IProbe (int, int, int&, MPI_Status&) {}

void Parallel::Comm_dup (MPI_Comm, MPI_Comm&) {}

void Parallel::ReduceDoubleMax (double&) {}
void Parallel::ReduceDoubleMin (double&) {}
void Parallel::ReduceDoubleSum (double&) {}

void Parallel::ReduceDoubleMax (double&,int) {}
void Parallel::ReduceDoubleMin (double&,int) {}
void Parallel::ReduceDoubleSum (double&,int) {}

void Parallel::ReduceDoubleMax (double*,int) {}
void Parallel::ReduceDoubleMin (double*,int) {}
void Parallel::ReduceDoubleSum (double*,int) {}

void Parallel::ReduceDoubleMax (double*,int,int) {}
void Parallel::ReduceDoubleMin (double*,int,int) {}
void Parallel::ReduceDoubleSum (double*,int,int) {}

void Parallel::ReduceLongAnd (long&) {}
void Parallel::ReduceLongSum (long&) {}
void Parallel::ReduceLongMax (long&) {}
void Parallel::ReduceLongMin (long&) {}

void Parallel::ReduceLongAnd (long&,int) {}
void Parallel::ReduceLongSum (long&,int) {}
void Parallel::ReduceLongMax (long&,int) {}
void Parallel::ReduceLongMin (long&,int) {}

void Parallel::ReduceLongAnd (long*,int) {}
void Parallel::ReduceLongSum (long*,int) {}
void Parallel::ReduceLongMax (long*,int) {}
void Parallel::ReduceLongMin (long*,int) {}

void Parallel::ReduceLongAnd (long*,int,int) {}
void Parallel::ReduceLongSum (long*,int,int) {}
void Parallel::ReduceLongMax (long*,int,int) {}
void Parallel::ReduceLongMin (long*,int,int) {}

void Parallel::ReduceIntSum (int&) {}
void Parallel::ReduceIntMax (int&) {}
void Parallel::ReduceIntMin (int&) {}

void Parallel::ReduceIntSum (int&,int) {}
void Parallel::ReduceIntMax (int&,int) {}
void Parallel::ReduceIntMin (int&,int) {}

void Parallel::ReduceIntSum (int*,int) {}
void Parallel::ReduceIntMax (int*,int) {}
void Parallel::ReduceIntMin (int*,int) {}

void Parallel::ReduceIntSum (int*,int,int) {}
void Parallel::ReduceIntMax (int*,int,int) {}
void Parallel::ReduceIntMin (int*,int,int) {}

void Parallel::ReduceBoolAnd (bool&) {}
void Parallel::ReduceBoolOr  (bool&) {}

void Parallel::ReduceBoolAnd (bool&,int) {}
void Parallel::ReduceBoolOr  (bool&,int) {}
//
// Here so we don't need to include <Utility.H> in <Parallel.H>.
//


//from Utility.cpp in BoxLib TODO proper credit
static
double
get_initial_wall_clock_time ()
{
    struct timeval tp;

    if (gettimeofday(&tp, 0) != 0)
    {
        std::cout << "get_time_of_day(): gettimeofday() failed" << std::endl;
        std::abort();
    }

    return tp.tv_sec + tp.tv_usec/1000000.0;
}

//
// Attempt to guarantee wsecond() gets initialized on program startup.
//
double Initial_Wall_Clock_Time = get_initial_wall_clock_time();

double
wsecond (double* t=0)
{
    struct timeval tp;

    gettimeofday(&tp,0);

    double dt = tp.tv_sec + tp.tv_usec/1000000.0 - Initial_Wall_Clock_Time;

    if (t != 0)
        *t = dt;

    return dt;
}
//END


double
Parallel::second ()
{
  return wsecond();
}

#endif
//
// This function is the same whether or not we're using MPI.
//
int
Parallel::SeqNum ()
{
    const int BEG = 1000;
    const int END = 9000;

    static int seqno = BEG;

    int result = seqno++;

    if (seqno > END) seqno = BEG;

    return result;
}



void
Parallel::ReadAndBcastFile (const std::string& filename,
                                      std::vector<char>&       charBuf)
{
    enum { IO_Buffer_Size = 40960 * 32 };


    typedef char Setbuf_Char_Type;


    std::vector<Setbuf_Char_Type> io_buffer(IO_Buffer_Size);

    int fileLength = 0, fileLengthPadded;

    std::ifstream iss;

    if (Parallel::IOProcessor())
    {
        iss.rdbuf()->pubsetbuf(io_buffer.data(), io_buffer.size());
        iss.open(filename.c_str(), std::ios::in);
        if (!iss.good())
        {
            std::abort();
        }
        iss.seekg(0, std::ios::end);
        fileLength = iss.tellg();
        iss.seekg(0, std::ios::beg);
    }
    Parallel::Bcast(&fileLength, 1,
                              Parallel::IOProcessorNumber());
    fileLengthPadded = fileLength + 1;
    fileLengthPadded += fileLengthPadded % 8;
    charBuf.resize(fileLengthPadded);
    if (Parallel::IOProcessor())
    {
        iss.read(charBuf.data(), fileLength);
        iss.close();
    }
    Parallel::Bcast(charBuf.data(), fileLengthPadded,
                              Parallel::IOProcessorNumber());
    charBuf[fileLength] = '\0';
}

