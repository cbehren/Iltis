#ifndef PARALLEL_H
#define PARALLEL_H

//This is a stripped/modified version of the ParallelDescriptor module as included in the BoxLib library, see https://github.com/BoxLib-Codes/BoxLib 
//For details of the license, see license.txt

//this is a way of making the MPI stuff more abstract. We rarely need to call MPI routines, because we can use the wrappers provided by this class. 


#include "mpi-wrapper.H"
#include <vector>
#include <string>

namespace Parallel
{
    
}

//
// Functions used for implementing parallelism.
//

namespace Parallel
{
    class Message
    {
    public:

	Message () :
            m_finished(true),
            m_type(MPI_DATATYPE_NULL),
            m_req(MPI_REQUEST_NULL) {}
	Message (MPI_Request req_, MPI_Datatype type_) :
            m_finished(false),
            m_type(type_),
            m_req(req_) {}
	Message (MPI_Status stat_, MPI_Datatype type_) :
            m_finished(true),
            m_type(type_),
            m_req(MPI_REQUEST_NULL), m_stat(stat_) {}
	void wait ();
	bool test ();
	size_t count () const;
	int tag () const;
	int pid () const;
	MPI_Datatype type () const { return m_type; }
	MPI_Request  req () const { return m_req; }

    private:

	bool               m_finished;
	MPI_Datatype       m_type;
	MPI_Request        m_req;
	mutable MPI_Status m_stat;
    };
    //
    // Perform any needed parallel initialization.  This MUST be the
    // first routine in this class called from within a program.
    //
    void StartParallel (int*    argc = 0,
			char*** argv = 0,
                        MPI_Comm mpi_comm = MPI_COMM_WORLD,
                        int required = 0
                       );
    //
    // Perform any needed parallel finalization.  This MUST be the
    // last routine in this class called from within a program.
    //
    void EndParallel ();
    //
    // Returns processor number of calling program.
    //
    extern int m_MyId;
    inline int
    MyProc ()
    {
        return m_MyId;
    }
    //
    // Returns number of CPUs involved in the computation.
    //
    extern int m_nProcs;
    inline int
    NProcs ()
    {
        return m_nProcs;
    }
    //
    // The CPU number of the I/O Processor.
    //
    extern const int ioProcessor;
    inline int
    IOProcessorNumber ()
    {
        return ioProcessor;
    }
    //
    // Is this CPU the I/O Processor?
    //
    inline bool
    IOProcessor ()
    {
         return MyProc() == IOProcessorNumber();
    }
    //
    //  Parallel Communicator, probably MPI_COMM_WORLD
    //
    extern MPI_Comm m_comm;
    inline MPI_Comm Communicator ()
    {
        return m_comm;
    }

    void Barrier ();
    void Barrier (MPI_Comm comm);

    void Test (MPI_Request& request, int& flag, MPI_Status& status);

    void Comm_dup (MPI_Comm comm, MPI_Comm& newcomm);
    //
    // Issue architecture specific Abort.
    //
    void Abort ();
    //
    // Abort with specified error code.
    //
    void Abort (int errorcode);
    //
    // ErrorString return string associated with error internal error condition
    //
    const char* ErrorString (int errcode);
    //
    // Returns wall-clock seconds since start of execution.
    //
    double second ();
    //
    // And-wise boolean reduction.
    //
    void ReduceBoolAnd (bool& rvar);
    //
    // And-wise boolean reduction to specified cpu.
    //
    void ReduceBoolAnd (bool& rvar, int cpu);
    //
    // Or-wise boolean reduction.
    //
    void ReduceBoolOr  (bool& rvar);
    //
    // Or-wise boolean reduction to specified cpu.
    //
    void ReduceBoolOr  (bool& rvar, int cpu);
    //
    // double sum reduction.
    //
    void ReduceDoubleSum (double& rvar);

    void ReduceDoubleSum (double* rvar, int cnt);
    //
    // double sum reduction to specified cpu.
    //
    void ReduceDoubleSum (double& rvar, int cpu);

    void ReduceDoubleSum (double* rvar, int cnt, int cpu);
    //
    // double max reduction.
    //
    void ReduceDoubleMax (double& rvar);

    void ReduceDoubleMax (double* rvar, int cnt);
    //
    // double max reduction to specified cpu.
    //
    void ReduceDoubleMax (double& rvar, int cpu);

    void ReduceDoubleMax (double* rvar, int cnt, int cpu);
    //
    // double min reduction.
    //
    void ReduceDoubleMin (double& rvar);

    void ReduceDoubleMin (double* rvar, int cnt);
    //
    // double min reduction to specified cpu.
    //
    void ReduceDoubleMin (double& rvar, int cpu);

    void ReduceDoubleMin (double* rvar, int cnt, int cpu);
    //
    // Integer sum reduction.
    //
    void ReduceIntSum (int& rvar);

    void ReduceIntSum (int* rvar, int cnt);
    //
    // Integer sum reduction to specified cpu.
    //
    void ReduceIntSum (int& rvar, int cpu);

    void ReduceIntSum (int* rvar, int cnt, int cpu);
    //
    // Integer max reduction.
    //
    void ReduceIntMax (int& rvar);

    void ReduceIntMax (int* rvar, int cnt);
    //
    // Integer max reduction to specified cpu.
    //
    void ReduceIntMax (int& rvar, int cpu);

    void ReduceIntMax (int* rvar, int cnt, int cpu);
    //
    // Integer min reduction.
    //
    void ReduceIntMin (int& rvar);

    void ReduceIntMin (int* rvar, int cnt);
    //
    // Integer min reduction to specified cpu.
    //
    void ReduceIntMin (int& rvar, int cpu);

    void ReduceIntMin (int* rvar, int cnt, int cpu);
    //
    // Long sum reduction.
    //
    void ReduceLongSum (long& rvar);

    void ReduceLongSum (long* rvar, int cnt);
    //
    // Long sum reduction to specified cpu.
    //
    void ReduceLongSum (long& rvar, int cpu);

    void ReduceLongSum (long* rvar, int cnt, int cpu);
    //
    // Long max reduction.
    //
    void ReduceLongMax (long& rvar);

    void ReduceLongMax (long* rvar, int cnt);
    //
    // Long max reduction to specified cpu.
    //
    void ReduceLongMax (long& rvar, int cpu);

    void ReduceLongMax (long* rvar, int cnt, int cpu);
    //
    // Long min reduction.
    //
    void ReduceLongMin (long& rvar);

    void ReduceLongMin (long* rvar, int cnt);
    //
    // Long min reduction to specified cpu.
    //
    void ReduceLongMin (long& rvar, int cpu);

    void ReduceLongMin (long* rvar, int cnt, int cpu);
    //
    // Long and-wise reduction.
    //
    void ReduceLongAnd (long& rvar);

    void ReduceLongAnd (long* rvar, int cnt);
    //
    // Long and-wise reduction to specified cpu.
    //
    void ReduceLongAnd (long& rvar, int cpu);

    void ReduceLongAnd (long* rvar, int cnt, int cpu);
    //
    // Parallel gather.
    //
    void Gather (double* sendbuf,
                 int   sendcount,
                 double* recvbuf,
                 int   root);
    void Gatherv (double* sendbuf,
                  int sendcount,
                  int*   recvcounts,
                  double* recvbuf,
                  int   root);
    void Gatherv (int* sendbuf,
                  int sendcount,
                  int*   recvcounts,
                  int* recvbuf,
                  int   root);
    void Gatherv (long* sendbuf,
                  int sendcount,
                  int*   recvcounts,
                  long* recvbuf,
                  int   root);
    //
    // Returns sequential message sequence numbers in range 1000-9000.
    //
    int SeqNum ();

    template <class T> Message Asend(const T*, size_t n, int pid, int tag);
    template <class T> Message Asend(const T*, size_t n, int pid, int tag, MPI_Comm comm);
    template <class T> Message Asend(const std::vector<T>& buf, int pid, int tag);

    template <class T> Message Arecv(T*, size_t n, int pid, int tag);
    template <class T> Message Arecv(T*, size_t n, int pid, int tag, MPI_Comm comm);
    template <class T> Message Arecv(std::vector<T>& buf, int pid, int tag);

    template <class T> Message Send(const T* buf, size_t n, int dst_pid, int tag);
    template <class T> Message Send(const std::vector<T>& buf, int dst_pid, int tag);

    template <class T> Message Recv(T*, size_t n, int pid, int tag);
    template <class T> Message Recv(std::vector<T>& t, int pid, int tag);

    template <class T> void Bcast(T*, size_t n, int root = 0);

    template <class Op, class T> T Reduce(const T& t);

    template <class T, class T1> void Scatter(T*, size_t n, const T1*, size_t n1, int root);

    template <class T, class T1> void Gather(const T*, size_t n, T1*, size_t n1, int root);
    template <class T> std::vector<T> Gather(const T&, int root);


    void MPI_Error(const char* file, int line, const char* msg, int rc);

    void ReadAndBcastFile(const std::string &filename, std::vector<char> &charBuf);
    void IProbe(int src_pid, int tag, int &mflag, MPI_Status &status);
}

#define BL_MPI_REQUIRE(x)						\
do									\
{									\
  if ( int l_status_ = (x) )						\
    {									\
      Parallel::MPI_Error(__FILE__,__LINE__,#x, l_status_);   \
    }									\
}									\
while ( false )

#if USE_MPI
template <class T>
Parallel::Message
Parallel::Asend (const T* buf,
                           size_t   n,
                           int      dst_pid,
                           int      tag)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(buf),
                              n,
                              Mpi_typemap<T>::type(),
                              dst_pid,
                              tag,
                              Communicator(),
                              &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
Parallel::Message
Parallel::Asend (const T* buf,
                           size_t   n,
                           int      dst_pid,
                           int      tag,
                           MPI_Comm comm)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(buf),
                              n,
                              Mpi_typemap<T>::type(),
                              dst_pid,
                              tag,
                              comm,
                              &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
Parallel::Message
Parallel::Asend (const std::vector<T>& buf,
                           int                   dst_pid,
                           int                   tag)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(&buf[0]),
                              buf.size(),
                              Mpi_typemap<T>::type(),
                              dst_pid,
                              tag,
                              Communicator(),
                              &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
Parallel::Message
Parallel::Send (const T* buf,
                          size_t   n,
                          int      dst_pid,
                          int      tag)
{
    BL_MPI_REQUIRE( MPI_Send(const_cast<T*>(buf),
                             n,
                             Mpi_typemap<T>::type(),
                             dst_pid,
                             tag,
                             Communicator()) );
    return Message();
}

template <class T>
Parallel::Message
Parallel::Send (const std::vector<T>& buf,
                          int                   dst_pid,
                          int                   tag)
{
    BL_MPI_REQUIRE( MPI_Send(const_cast<T*>(&buf[0]),
                             buf.size(),
                             Mpi_typemap<T>::type(),
                             dst_pid,
                             tag,
                             Communicator()) );
    return Message();
}

template <class T>
Parallel::Message
Parallel::Arecv (T*       buf,
                           size_t   n,
                           int      src_pid,
                           int      tag)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(buf,
                              n,
                              Mpi_typemap<T>::type(),
                              src_pid,
                              tag,
                              Communicator(),
                              &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
Parallel::Message
Parallel::Arecv (T*       buf,
                           size_t   n,
                           int      src_pid,
                           int      tag,
                           MPI_Comm comm)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(buf,
                              n,
                              Mpi_typemap<T>::type(),
                              src_pid,
                              tag,
                              comm,
                              &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
Parallel::Message
Parallel::Arecv (std::vector<T>& buf,
                           int             src_pid,
                           int             tag)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(&buf[0],
                              buf.size(),
                              Mpi_typemap<T>::type(),
                              src_pid,
                              tag,
                              Communicator(),
                              &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
Parallel::Message
Parallel::Recv (T*     buf,
                          size_t n,
                          int    src_pid,
                          int    tag)
{
    MPI_Status stat;
    BL_MPI_REQUIRE( MPI_Recv(buf,
                             n,
                             Mpi_typemap<T>::type(),
                             src_pid,
                             tag,
                             Communicator(),
                             &stat) );
    return Message(stat, Mpi_typemap<T>::type());
}

template <class T>
Parallel::Message
Parallel::Recv (std::vector<T>& buf,
                          int             src_pid,
                          int             tag)
{
    MPI_Status stat;
    BL_MPI_REQUIRE( MPI_Recv(&buf[0],
                             buf.size(),
                             Mpi_typemap<T>::type(),
                             src_pid,
                             tag,
                             Communicator(),
                             &stat) );
    return Message(stat, Mpi_typemap<T>::type());
}

template <class Op, class T>
T
Parallel::Reduce (const T& t)
{
    T recv;
    BL_MPI_REQUIRE( MPI_Allreduce(const_cast<T*>(&t),
                                  &recv,
                                  1,
                                  Mpi_typemap<T>::type(),
                                  Op::op(),
                                  Communicator()) );
    return recv;
}

template <class T>
void
Parallel::Bcast (T*     t,
                           size_t n,
                           int    root)
{
    BL_MPI_REQUIRE( MPI_Bcast(t,
                              n,
                              Mpi_typemap<T>::type(),
                              root,
                              Communicator()) );
}

template <class T, class T1>
void
Parallel::Gather (const T* t,
                            size_t   n,
                            T1*      t1,
                            size_t   n1,
                            int      root)
{
    BL_MPI_REQUIRE( MPI_Gather(const_cast<T*>(t),
                               n,
                               Mpi_typemap<T>::type(),
                               t1,
                               n1,
                               Mpi_typemap<T1>::type(),
                               root,
                               Communicator()) );
}

template <class T>
std::vector<T>
Parallel::Gather (const T& t, int root)
{
    std::vector<T> resl;
    if ( root == MyProc() ) resl.resize(NProcs());
    BL_MPI_REQUIRE( MPI_Gather(const_cast<T*>(&t),
                               1,
                               Mpi_typemap<T>::type(),
                               &resl[0],
                               1,
                               Mpi_typemap<T>::type(),
                               root,
                               Communicator()) );
    return resl;
}

template <class T, class T1>
void
Parallel::Scatter (T*        t,
                             size_t    n,
                             const T1* t1,
                             size_t    n1,
                             int       root)
{
    BL_MPI_REQUIRE( MPI_Scatter(const_cast<T1*>(t1),
                                n1,
                                Mpi_typemap<T1>::type(),
                                t,
                                n,
                                Mpi_typemap<T>::type(),
                                root,
                                Communicator()) );
}

#else

namespace Parallel
{
template <class T>
Message
Asend(const T* buf, size_t n, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Asend(const T* buf, size_t n, int dst_pid, int tag, MPI_Comm comm)
{
    return Message();
}

template <class T>
Message
Asend(const std::vector<T>& buf, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Send(const T* buf, size_t n, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Send(const std::vector<T>& buf, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Arecv(T* buf, size_t n, int src_pid, int tag)
{
    return Message();
}

template <class T>
Message
Arecv(T* buf, size_t n, int src_pid, int tag, MPI_Comm comm)
{
    return Message();
}

template <class T>
Message
Arecv(std::vector<T>& buf, int src_pid, int tag)
{
    return Message();
}

template <class T>
Message
Recv(T* buf, size_t n, int src_pid, int tag)
{
    return Message();
}

template <class T>
Message
Recv(std::vector<T>& buf, int src_pid, int tag)
{
    return Message();
}

template <class Op, class T>
T
Reduce(const T& t)
{
    return t;
}

template <class T>
void
Bcast(T* t, size_t n, int root)
{}

template <class T, class T1>
void
Gather(const T* t, size_t n, T1* t1, size_t n1, int root)
{}

template <class T>
std::vector<T>
Gather(const T& t, int root)
{
    std::vector<T> resl(1);
    resl[0] = t;
    return resl;
}

template <class T, class T1>
void
Scatter(T* t, size_t n, const T1* t1, size_t n1, int root)
{}

}
#endif

#endif 
