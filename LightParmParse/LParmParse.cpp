
#include <typeinfo>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <vector>
#include <list>

#include "LParmParse.H"

//This is a stripped/modified version of ParmParse as included in the BoxLib library, see https://github.com/BoxLib-Codes/BoxLib 
//For details of the license, see license.txt

//
// Used by constructor to build table.
//
LParmParse::PP_entry::PP_entry (const std::string& name,
                    const std::list<std::string>& vals)
    :
    m_name(name),
    m_table(0),
    m_queried(false)
{
    for ( std::list<std::string>::const_iterator li = vals.begin(), End = vals.end(); li != End; ++li )
    {
	m_vals.push_back(*li);
    }
}

LParmParse::PP_entry::PP_entry (const std::string& name,
                    const std::string&            val)
    :
    m_name(name),
    m_table(0),
    m_queried(false)
{
    m_vals.push_back(val);
}

LParmParse::PP_entry::PP_entry (const std::string& name,
                    const std::list<PP_entry>&   table)
    :
    m_name(name),
    m_table(new Table(table)),
    m_queried(false)
{
}

LParmParse::PP_entry::PP_entry (const PP_entry& pe)
    : m_name(pe.m_name),
      m_vals(pe.m_vals),
      m_table(0),
      m_queried(pe.m_queried)
{
    if ( pe.m_table )
    {
	m_table = new Table(*pe.m_table);
    }
}

LParmParse::PP_entry::~PP_entry ()
{
    delete m_table;
}
    
LParmParse::PP_entry&
LParmParse::PP_entry::operator= (const PP_entry& pe)
{
    if ( &pe == this ) return *this;
    m_name = pe.m_name;
    m_vals = pe.m_vals;
    m_table = 0;
    m_queried = pe.m_queried;
    if ( pe.m_table )
    {
	m_table = new Table(*pe.m_table);
    }
    return *this;
}

std::string
LParmParse::PP_entry::print () const {
    std::stringstream t;
    t << m_name << " = ";
    int n = m_vals.size();
    for ( int i = 0; i < n; i++)
    {
	t << m_vals[i];
	if ( i < n-1 ) t << " ";
    }
    return t.str();
}

std::ostream&
operator<< (std::ostream& os, const LParmParse::PP_entry& pp)
{
    os << pp.m_name << "(nvals = " << pp.m_vals.size() << ") " << " :: [";
    int n = pp.m_vals.size();
    for ( int i = 0; i < n; i++ )
    {
	os << pp.m_vals[i];
	if ( i < n-1 ) os << ", ";
    }
    os << "]";

    if ( !os )
    {
        std::abort();
    }
    return os;
}

namespace
{
enum PType
{
    pDefn,
    pValue,
    pEQ_sign,
    pOpenBracket,
    pCloseBracket,
    pEOF
};

template <class T>
bool
is (const std::string& str, T& val)
{
    std::istringstream s(str);
    s >> val;
    if ( !s ) return false;
    return true;
}

template <>
bool
is (const std::string& str, std::string& val)
{
    val = str;
    return true;
}

template <> 
bool
is (const std::string& str, bool& val)
{
    if ( str == "true" || str == "t" )
    {
        val = true;
        return true;
    }
    if ( str == "false" || str == "f" )
    {
        val = false;
        return true;
    }
    int int_val;
    if ( is(str, int_val) )
    {
        val = int_val != 0;
        return true;
    }
    double dbl_val;
    if ( is(str, dbl_val) )
    {
        val = dbl_val != 0;
        return true;
    }
    return false;
}

LParmParse::Table g_table;
typedef std::list<LParmParse::PP_entry>::iterator list_iterator;
typedef std::list<LParmParse::PP_entry>::const_iterator const_list_iterator;

template <class T> const char* tok_name(const T&) { return typeid(T).name(); }
template <class T> const char* tok_name(std::vector<T>&) { return tok_name(T());}

//
// Simple lexical analyser.
//

enum lexState
{
    START,
    STRING,
    QUOTED_STRING,
    IDENTIFIER,
    LIST
};

const char* const
state_name[] =
{
   "START",
   "STRING",
   "QUOTED_STRING",
   "IDENTIFIER",
   "LIST"
};

void
eat_garbage (const char*& str)
{
    for (;;)
    {
	if ( *str == 0 ) break;
        else if ( *str == '#' )
        {
            while ( *str && *str != '\n' )
	    {
		str++;
	    }
	    continue;
        }
        else if ( isspace(*str) )
	{
	    str++;
	}
        else
	{
            break;
	}
    }
}

PType
getToken (const char*& str,
	  std::string& ostr)
{
   //
   // Eat white space and comments.
   //
   eat_garbage(str);
   //
   // Check for end of file.
   //
   if ( *str == 0 )
   {
       return pEOF;
   }
   //
   // Start token scan.
   //
   lexState state = START;
   int      pcnt  = 0; // Tracks nested parens
   while (true)
   {
       char ch = *str;
       if ( ch == 0 )
       {
	   std::abort();
       }
       switch (state)
       {
       case START:
           if ( ch == '=' )
           {
               ostr += ch; str++;
               return pEQ_sign;
           }
           else if ( ch == '"' )
           {
               str++;
               state = QUOTED_STRING;
           }
	   else if ( ch == '(' )
	   {
	       ostr += ch; str++; pcnt = 1;
	       state = LIST;
	   }
	   else if ( ch == '{' )
	   {
	       str++;
	       return pOpenBracket;
	   }
	   else if ( ch == '}' )
	   {
	       str++;
	       return pCloseBracket;
	   }
           else if ( isalpha(ch) )
           {
               ostr += ch; str++;
               state = IDENTIFIER;
           }
           else
           {
               ostr += ch; str++;
               state = STRING;
           }
           break;
       case IDENTIFIER:
           if ( isalnum(ch) || ch == '_' || ch == '.' || ch == '[' || ch == ']' || ch == '+' || ch == '-' )
           {
               ostr += ch; str++;
           }
           else if ( isspace(ch) || ch == '=' )
           {
               return pDefn;
           }
           else
           {
               ostr += ch; str++;
               state = STRING;
           }
           break;
       case LIST:
	   if ( ch == '(' )
	   {
	       ostr += ch; str++; pcnt++;
	   }
	   else if ( ch == ')' )
	   {
	       ostr += ch; str++; pcnt--;
	       if ( pcnt == 0 )
               {
		   return pValue;
               }
	   }
	   else
	   {
	       ostr += ch; str++;
	   }
	   break;
       case STRING:
           if ( isspace(ch) || ch == '=' )
           {
               return pValue;
           }
           else
           {
               ostr += ch; str++;
           }
           break;
       case QUOTED_STRING:
           if ( ch == '"' )
           {
               str++;
               return pValue;
           }
           else
           {
               ostr += ch; str++;
           }
           break;
       default:
           std::abort();
            }
   }

}


//
// Keyword aware string comparison.
//

static
bool
ppfound (const std::string& keyword,
	 const LParmParse::PP_entry& pe,
	 bool recordQ)
{
    return (recordQ == (pe.m_table!=0)) && (keyword == pe.m_name);
}

//
// Return the index of the n'th occurence of a parameter name,
// except if n==-1, return the index of the last occurence.
// Return 0 if the specified occurence does not exist.
//

const LParmParse::PP_entry*
ppindex (const LParmParse::Table& table,
	 int         n,
	 const std::string& name,
	 bool recordQ)
{
    const LParmParse::PP_entry* fnd = 0;

    if ( n == LParmParse::LAST )
    {
        //
        // Search from back of list.
        //
        for (std::list<LParmParse::PP_entry>::const_reverse_iterator li = table.rbegin(), REnd = table.rend(); li != REnd; ++li)
        {
            if ( ppfound(name, *li, recordQ) )
            {
                fnd = &*li;
                break;
            }
        }
    }
    else
    {
        for ( const_list_iterator li =table.begin(), End = table.end(); li != End; ++li )
        {
            if ( ppfound(name, *li, recordQ) )
            {
                fnd = &*li;
                if ( --n < 0 )
		{
                    break;
		}
            }
        }
        if ( n >= 0)
	{
            fnd = 0;
	}
    }

    if ( fnd )
    {
        //
        // Found an entry; mark all occurences of name as used.
        //
        for ( const_list_iterator li = table.begin(), End = table.end(); li != End; ++li )
	{
            if ( ppfound(name, *li, recordQ) )
	    {
                li->m_queried = true;
	    }
	}
    }
    return fnd;
}

void
bldTable (const char*& str, std::list<LParmParse::PP_entry>& tab);

void ReadFile(const std::string& filename,std::vector<char>& charBuf)
{
    int fileLength = 0, fileLengthPadded;
    std::ifstream iss;
    iss.open(filename.c_str(), std::ios::in);
    if (!iss.good())
        {
            std::abort();
        }
    iss.seekg(0, std::ios::end);
    fileLength = iss.tellg();
    iss.seekg(0, std::ios::beg);
    fileLengthPadded = fileLength + 1;
    fileLengthPadded += fileLengthPadded % 8;
    charBuf.resize(fileLengthPadded);
    iss.read(charBuf.data(), fileLength);
    iss.close();
    charBuf[fileLength] = '\0';
}

static void
read_file (const char*                     fname,
	   std::list<LParmParse::PP_entry>& tab)
{
    //
    // Space for input file if it exists.
    //
    if ( fname != 0 && fname[0] != 0 )
    {
	std::vector<char> fileCharPtr;
	std::string filename = fname;
    ReadFile(fname,fileCharPtr);
	const char* b = fileCharPtr.data();
        bldTable(b, tab);
    }
}

void
addDefn (std::string&         def,
	 std::list<std::string>&   val,
	 std::list<LParmParse::PP_entry>& tab)
{
    static const std::string FileKeyword("FILE");
    //
    // Check that defn exists.
    //
    if ( def.empty() )
    {
        val.clear();
        return;
    }
    //
    // Check that it has values.
    //
    if ( val.empty() )
    {
        std::cerr << "LParmParse::addDefn(): no values for definition " << def << "\n";
        std::abort();
    }
    //
    // Check if this defn is a file include directive.
    //
    if ( def == FileKeyword && val.size() == 1 )
    {
        //
        // Read file and add to this table.
        //
        const char* fname = val.front().c_str();
        read_file(fname, tab);
    }
    else
    {
        tab.push_back(LParmParse::PP_entry(def,val));
    }
    val.clear();
    def = std::string();
}

void
addTable (std::string& def,
	  LParmParse::Table& val,
	  std::list<LParmParse::PP_entry>& tab)
{
    if ( def.empty() )
    {
        val.clear();
        return;
    }
    //
    // Check that it has values.
    //
    if ( val.empty() )
    {
        std::cerr << "LParmParse::addTable(): no values for Table " << def << "\n";
	std::abort();
    }
    tab.push_back(LParmParse::PP_entry(def, val));
    val.clear();
    def = std::string();
}

void
bldTable (const char*&                    str,
	  std::list<LParmParse::PP_entry>& tab)
{
    std::string            cur_name;
    std::list<std::string> cur_list;
    LParmParse::Table       cur_table;
    std::string            tmp_str;

    for (;;)
    {
        std::string tokname;

	PType token = getToken(str,tokname);

	switch (token)
	{
	case pCloseBracket:
	    if ( !cur_name.empty() && cur_list.empty() )
	    {
		std::abort();
	    }
	case pEOF:
	    addDefn(cur_name,cur_list,tab);
	    return;
	case pOpenBracket:
	    if ( cur_name.empty() )
	    {
		std::abort();
	    }
	    if ( !cur_list.empty() )
	    {
		tmp_str = cur_list.back();
		cur_list.pop_back();
		addDefn(cur_name, cur_list, tab);
		cur_name = tmp_str;
	    }
	    bldTable(str, cur_table);
	    addTable(cur_name, cur_table, tab);
	    break;
	case pEQ_sign:
	    if ( cur_name.empty() )
	    {
		std::abort();
	    }
	    if ( !cur_list.empty() )
	    {
		tmp_str = cur_list.back();
		cur_list.pop_back();
		addDefn(cur_name,cur_list,tab);
		cur_name = tmp_str;
	    }
	    //
	    // Read one too far, remove last name on list.
	    //
	    break;
	case pDefn:
	    if ( cur_name.empty() )
	    {
		cur_name = tokname;
		break;
	    }
	    //
	    // Otherwise, fall through, this may be a string.
	    //
	case pValue:
	    if ( cur_name.empty() )
	    {
		std::string msg("LParmParse::bldTable(): value with no defn: ");
		msg += tokname;
		std::abort();
	    }
	    cur_list.push_back(tokname);
	    break;
	}
    }
}

namespace
{
template <class T>
bool
squeryval (const LParmParse::Table& table,
	   const std::string& name,
	   T&           ptr,
	   int          ival,
	   int          occurence)
{
    //
    // Get last occurrance of name in table.
    //
    const LParmParse::PP_entry* def = ppindex(table, occurence, name, false);
    if ( def == 0 )
    {
        return false;
    }
    //
    // Does it have ival values?
    //
    if ( ival >= (int)def->m_vals.size() )
    {
        std::cerr << "LParmParse::queryval no value number"
                  << ival << " for ";
        if ( occurence ==  LParmParse::LAST )
	{
            std::cerr << "last occurence of ";
	}
        else
	{
            std::cerr << " occurence " << occurence << " of ";
	}
        std::cerr << def->m_name << '\n' << *def << '\n';
        std::abort();
    }

    const std::string& valname = def->m_vals[ival];

    bool ok = is(valname, ptr);
    if ( !ok )
    {
        std::cerr << "LParmParse::queryval type mismatch on value number "
                  << ival << " of " << '\n';
        if ( occurence == LParmParse::LAST )
	{
            std::cerr << " last occurence of ";
	}
        else
	{
            std::cerr << " occurence number " << occurence << " of ";
	}
        std::cerr << def->m_name << '\n';
        std::cerr << " Expected an \""
                  << tok_name(ptr)
                  << "\" type  which can't be parsed from the string \""
                  << valname << "\"\n"
		  << *def << '\n';
        std::abort();
    }
    return true;
}

template <class T>
void
sgetval (const LParmParse::Table& table,
	 const std::string& name,
	 T&           ptr,
	 int          ival,
	 int          occurence)
{
    if ( squeryval(table, name,ptr,ival,occurence) == 0 )
    {
        std::cerr << "LParmParse::getval ";
        if ( occurence >= 0 )
	{
            std::cerr << "occurence number "
                      << occurence
                      << " of ";
	}

        std::cerr << "LParmParse::getval(): "
                  << name
                  << " not found in table"
                  << '\n';
        LParmParse::dumpTable(std::cerr);
        std::abort();
    }
}

template <class T>
bool
squeryarr (const LParmParse::Table& table,
	   const std::string& name,
	   std::vector<T>&    ptr,
	   int          start_ix,
	   int          num_val,
	   int          occurence)
{
    //
    // Get last occurrance of name in table.
    //
    const LParmParse::PP_entry *def = ppindex(table,occurence, name, false);
    if ( def == 0 )
    {
        return false;
    }
    //
    // Does it have sufficient number of values and are they all
    // the same type?
    //
    if ( num_val == LParmParse::ALL )
    {
	num_val = def->m_vals.size();
    }

    if ( num_val == 0 ) return true;

    int stop_ix = start_ix + num_val - 1;
    if ( (int)ptr.size() <= stop_ix )
    {
        ptr.resize(stop_ix + 1);
    }
    if ( stop_ix >= (int)def->m_vals.size() )
    {
        std::cerr << "LParmParse::queryarr too many values requested for";
        if ( occurence == LParmParse::LAST )
	{
            std::cerr << " last occurence of ";
	}
        else
	{
            std::cerr << " occurence " << occurence << " of ";
	}
        std::cerr << def->m_name << '\n' << *def << '\n';
        std::abort();
    }
    for ( int n = start_ix; n <= stop_ix; n++ )
    {
	const std::string& valname = def->m_vals[n];
	bool ok = is(valname, ptr[n]);
	if ( !ok )
	{
	    std::cerr << "LParmParse::queryarr type mismatch on value number "
		      <<  n << " of ";
	    if ( occurence == LParmParse::LAST )
	    {
		std::cerr << " last occurence of ";
	    }
	    else
	    {
		std::cerr << " occurence number " << occurence << " of ";
	    }
	    std::cerr << def->m_name << '\n';
	    std::cerr << " Expected an \""
		      << tok_name(ptr)
		      << "\" type which can't be parsed from the string \""
		      << valname << "\"\n"
		      << *def << '\n';
	    std::abort();
	}
    }
    return true;
}

template <class T>
void
sgetarr (const LParmParse::Table& table,
	 const std::string&  name,
	 std::vector<T>&           ptr,
	 int          start_ix,
	 int          num_val,
	 int          occurence)
{
    if ( squeryarr(table,name,ptr,start_ix,num_val,occurence) == 0 )
    {
        std::cerr << "LParmParse::sgetarr ";
        if ( occurence >= 0 )
	{
            std::cerr << "occurence number " << occurence << " of ";
	}
        std::cerr << "LParmParse::sgetarr(): "
                  << name
                  << " not found in table"
                  << '\n';
        LParmParse::dumpTable(std::cerr);
        std::abort();
    }
}

template <class T>
void
saddval (const std::string&      name,
	 const T&                ptr)
{
	std::stringstream val;
	val << ptr;
	LParmParse::PP_entry entry(name,val.str());
	entry.m_queried=true;
	g_table.push_back(entry);
}


template <class T>
void
saddarr (const std::string&      name,
	 const std::vector<T>&   ptr)
{
	std::list<std::string> arr;
	for(unsigned int i = 0; i < ptr.size(); i++) {
		std::stringstream val;
		val << ptr[i];
		arr.push_back(val.str());
	}
	LParmParse::PP_entry entry(name,arr);
	entry.m_queried=true;
	g_table.push_back(entry);
}

}

//
// Initialize LParmParse.
//

bool initialized = false;

void
ppinit (int argc, char** argv, const char* parfile, LParmParse::Table& table)
{
    if ( parfile != 0 )
    {
	read_file(parfile, table);
    }

    if ( argc > 0 )
    {
        std::string argstr;
        const char SPACE = ' ';
        for ( int i = 0; i < argc; i++ )
        {
            argstr += argv[i];
            argstr += SPACE;
        }
        std::list<LParmParse::PP_entry> arg_table;
	const char* b = argstr.c_str();
        bldTable(b, arg_table);
        //
        // Append arg_table to end of existing table.
        //
        g_table.splice(table.end(), arg_table);
    }
    initialized = true;
}

}  // End of unnamed namespace.

std::string
LParmParse::prefixedName (const std::string& str) const
{
    if ( str.empty() )
    {
	std::abort();
    }
    if ( !m_pstack.top().empty())
    {
        return m_pstack.top() + '.' + str;
    }
    return str;
}

void
LParmParse::pushPrefix (const std::string& str)
{
    std::string s(str);
    if ( !s.empty() )
    {
	if ( !m_pstack.top().empty() )
	{
	    s = m_pstack.top() + "." + s;
	}
	m_pstack.push(s);
    }
}

void
LParmParse::popPrefix ()
{
    if ( m_pstack.size() <= 1 )
    {
	std::abort();
    }
    m_pstack.pop();
}

std::string
LParmParse::getPrefix() const
{
    return m_pstack.top();
}

LParmParse::LParmParse (const std::string& prefix)
    :
    m_table(g_table)
{
    m_pstack.push(prefix);
}

LParmParse::LParmParse (const Table& table)
    : m_table(table)
{
    m_pstack.push("");
}

LParmParse::Frame::Frame (LParmParse& pp, const std::string& pfix)
    :
    m_pp(pp), m_np(0)
{
    push(pfix);
    
}

LParmParse::Frame::~Frame ()
{
    
    while ( m_np )
    {
	pop();
    }
    
}

void
LParmParse::Frame::push (const std::string& str)
{
    m_pp.pushPrefix(str);
    m_np++;
}

void
LParmParse::Frame::pop ()
{
    
    m_pp.popPrefix();
    m_np--;
}

std::string
LParmParse::Frame::getPrefix () const
{
    return m_pp.getPrefix();
}

void
LParmParse::appendTable(LParmParse::Table& tab)
{
  g_table.splice(g_table.end(), tab);
}

static
bool
unused_table_entries_q (const LParmParse::Table& table)
{
    for ( const_list_iterator li = table.begin(), End = table.end(); li != End; ++li )
    {
	if ( li->m_table )
	{
	    if ( !li->m_queried )
	    {
		return true;
	    }
	    else
	    {
		return unused_table_entries_q(*li->m_table);
	    }
	}
	else if ( !li->m_queried )
	{
	    return true;
	}
    }
    return false;
}

static
void
finalize_table (std::string pfx, const LParmParse::Table& table)
{
    for ( const_list_iterator li = table.begin(), End = table.end(); li != End; ++li )
    {
	if ( li->m_table )
	{
	    if ( !li->m_queried )
	    {
		std::cout << "Record " << li->m_name << std::endl;
	    }
	    else
	    {
		finalize_table(pfx + "::" + li->m_name, *li->m_table);
	    }
	}
	else if ( !li->m_queried )
	{
	    std::cout << pfx << "::" << *li << std::endl;
	}
    }
}

void
LParmParse::Initialize (int         argc,
		       char**      argv,
		       const char* parfile)
{
    if ( initialized )
    {
	std::abort();
    }
    ppinit(argc, argv, parfile, g_table);

}

void
LParmParse::Finalize ()
{
    if ( unused_table_entries_q(g_table))
    {
	std::cout << "Unused LParmParse Variables:\n";
	finalize_table("[TOP]", g_table);
	std::cout << "done.\n";
	//
	// First loop through and delete all queried entries.
	//
    }
    g_table.clear();

    initialized = false;
}

void
LParmParse::dumpTable (std::ostream& os, bool prettyPrint)
{
    for ( const_list_iterator li = g_table.begin(), End = g_table.end(); li != End; ++li )
    {
	if(prettyPrint && li->m_queried) {
	    os << li->print() << std::endl;
	}
	else
	    os << *li << std::endl;
    }
}

int
LParmParse::countval (const char* name,
                     int         n) const
{
    //
    // First find n'th occurance of name in table.
    //
    const PP_entry* def = ppindex(m_table, n, prefixedName(name), false);
    return def == 0 ? 0 : def->m_vals.size();
}

// BOOL
void
LParmParse::getkth (const char* name,
                   int         k,
                   bool&        ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
LParmParse::get (const char* name,
                bool&        ptr,
                int ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
LParmParse::querykth (const char* name,
                     int         k,
                     bool&        ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
LParmParse::query (const char* name,
                  bool&        ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
LParmParse::add (const char* name,
                const bool  val)
{
    saddval(prefixedName(name),val);
}

// INT
void
LParmParse::getkth (const char* name,
                   int         k,
                   int&        ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
LParmParse::get (const char* name,
                int&        ptr,
                int ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
LParmParse::querykth (const char* name,
                     int         k,
                     int&        ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
LParmParse::query (const char* name,
                  int&        ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
LParmParse::add (const char* name,
                const int  val)
{
    saddval(prefixedName(name),val);
}

void
LParmParse::getktharr (const char* name,
                      int         k,
                      std::vector<int>& ptr,
                      int         start_ix,
                      int         num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
LParmParse::getarr (const char* name,
                   std::vector<int>& ptr,
                   int         start_ix,
                   int         num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
LParmParse::queryktharr (const char* name,
                        int         k,
                        std::vector<int>& ptr,
                        int         start_ix,
                        int         num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

int
LParmParse::queryarr (const char* name,
                     std::vector<int>& ptr,
                     int         start_ix,
                     int         num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

void
LParmParse::addarr (const char* name,
                const std::vector<int>&  ptr)
{
    saddarr(prefixedName(name),ptr);
}


// LONG
void
LParmParse::getkth (const char* name,
                   int         k,
                   long&       ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
LParmParse::get (const char* name,
                long&       ptr,
                int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
LParmParse::querykth (const char* name,
                     int         k,
                     long&       ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
LParmParse::query (const char* name,
                  long&       ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
LParmParse::add (const char* name,
                const long  val)
{
    saddval(prefixedName(name),val);
}

void
LParmParse::getktharr (const char* name,
                      int         k,
                      std::vector<long>& ptr,
                      int         start_ix,
                      int         num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}
void 
LParmParse::get_line_of_sight(const char* name,line_of_sight &los, int k)
{
    std::vector<double> temp;
    getktharr(name,k,temp,LParmParse::FIRST,3);
    for(int j=0;j<3;j++)
    {
        los.x[j] = temp[j];
    }
    
}


void
LParmParse::getarr (const char* name,
                   std::vector<long>& ptr,
                   int         start_ix,
                   int         num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
LParmParse::queryktharr (const char* name,
                        int         k,
                        std::vector<long>& ptr,
                        int         start_ix,
                        int         num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

int
LParmParse::queryarr (const char* name,
                     std::vector<long>& ptr,
                     int         start_ix,
                     int         num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

void
LParmParse::addarr (const char* name,
                const std::vector<long>&  ptr)
{
    saddarr(prefixedName(name),ptr);
}




// FLOAT
void
LParmParse::getkth (const char* name,
                   int         k,
                   float&      ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
LParmParse::get (const char* name,
                float&      ptr,
                int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
LParmParse::querykth (const char* name,
                     int         k,
                     float&      ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
LParmParse::query (const char* name,
                  float&      ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
LParmParse::add (const char* name,
                const float  val)
{
    saddval(prefixedName(name),val);
}

void
LParmParse::getktharr (const char*   name,
                      int           k,
                      std::vector<float>& ptr,
                      int           start_ix,
                      int           num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
LParmParse::getarr (const char*   name,
                   std::vector<float>& ptr,
                   int           start_ix,
                   int           num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
LParmParse::queryktharr (const char*   name,
                        int           k,
                        std::vector<float>& ptr,
                        int           start_ix,
                        int           num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix, num_val,k);
}

int
LParmParse::queryarr (const char*   name,
                     std::vector<float>& ptr,
                     int           start_ix,
                     int           num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

void
LParmParse::addarr (const char* name,
                const std::vector<float>&  ptr)
{
    saddarr(prefixedName(name),ptr);
}



// DOUBLE
void
LParmParse::getkth (const char* name,
                   int         k,
                   double&     ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
LParmParse::get (const char* name,
                double&     ptr,
                int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
LParmParse::querykth (const char* name,
                     int         k,
                     double&     ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
LParmParse::query (const char* name,
                  double&     ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
LParmParse::add (const char* name,
                const double  val)
{
    saddval(prefixedName(name),val);
}

void
LParmParse::getktharr (const char*    name,
                      int            k,
                      std::vector<double>& ptr,
                      int            start_ix,
                      int            num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
LParmParse::getarr (const char*    name,
                   std::vector<double>& ptr,
                   int            start_ix,
                   int            num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
LParmParse::queryktharr (const char*    name,
                        int            k,
                        std::vector<double>& ptr,
                        int            start_ix,
                        int            num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix, num_val,k);
}

int
LParmParse::queryarr (const char*    name,
                     std::vector<double>& ptr,
                     int            start_ix,
                     int            num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

void
LParmParse::addarr (const char* name,
                const std::vector<double>&  ptr)
{
    saddarr(prefixedName(name),ptr);
}



// STRING
void
LParmParse::getkth (const char* name,
                   int         k,
                   std::string&    ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
LParmParse::get (const char* name,
                std::string&    ptr,
                int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
LParmParse::querykth (const char* name,
                     int         k,
                     std::string&    ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
LParmParse::query (const char* name,
                  std::string&    ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
LParmParse::add (const char* name,
                const std::string&  val)
{
    saddval(prefixedName(name),val);
}

void
LParmParse::getktharr (const char*     name,
                      int             k,
                      std::vector<std::string>& ptr,
                      int             start_ix,
                      int             num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
LParmParse::getarr (const char*     name,
                   std::vector<std::string>& ptr,
                   int             start_ix,
                   int             num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
LParmParse::queryktharr (const char*     name,
                        int             k,
                        std::vector<std::string>& ptr,
                        int             start_ix,
                        int             num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix, num_val,k);
}

int
LParmParse::queryarr (const char*     name,
                     std::vector<std::string>& ptr,
                     int             start_ix,
                     int             num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

void
LParmParse::addarr (const char* name,
                const std::vector<std::string>&  ptr)
{
    saddarr(prefixedName(name),ptr);
}





//
// Return number of occurences of parameter name.
//

int
LParmParse::countname (const std::string& name) const
{
    int cnt = 0;
    for ( const_list_iterator li = m_table.begin(), End = m_table.end(); li != End; ++li )
    {
	if ( ppfound(prefixedName(name), *li, false) )
	{
	    cnt++;
	}
    }
    return cnt;
}

int
LParmParse::countRecords (const std::string& name) const
{
    int cnt = 0;
    for ( const_list_iterator li = m_table.begin(), End = m_table.end(); li != End; ++li )
    {
	if ( ppfound(prefixedName(name), *li, true) )
	{
	    cnt++;
	}
    }
    return cnt;
}

//
// Return true if name in table.
//

bool
LParmParse::contains (const char* name) const
{
    for ( const_list_iterator li = m_table.begin(), End = m_table.end(); li != End; ++li )
    {
       if ( ppfound(prefixedName(name), *li, false))
       {
           //
           // Found an entry; mark all occurences of name as used.
           //
           for ( const_list_iterator lli = m_table.begin(); lli != m_table.end(); ++lli )
	   {
               if ( ppfound(prefixedName(name), *lli, false) )
	       {
                   lli->m_queried = true;
	       }
	   }
           return true;
       }
    }
    return false;
}

LParmParse::Record
LParmParse::getRecord (const std::string& name, int n) const
{
    const PP_entry* pe = ppindex(m_table, n, prefixedName(name), true);
    if ( pe == 0 )
    {
	std::cerr << "LParmParse::getRecord: record " << name << " not found" << std::endl;
	std::abort();
    }
    return Record(LParmParse(*pe->m_table));
}


//
//
//

LParmParse::Record::Record ( const LParmParse& pp )
    : m_pp(pp)
{
}

const LParmParse*
LParmParse::Record::operator-> () const
{
    return &m_pp;
}

const LParmParse&
LParmParse::Record::operator* () const
{
    return m_pp;
}
