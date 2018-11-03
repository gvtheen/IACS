#include <stdio.h>
#include <limits.h>
#include "unistd.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "CIOMonitor.h"
#include "../Util/log.hpp"
#include "../GACatalyst.h"

using util::Log;

namespace CALCZJUT{

CIOMonitor::CIOMonitor(CParameter* m_parameter)
:m_pParameter(m_parameter)
{
    char abs_path_buff[PATH_MAX];
    if(NULL==getcwd(abs_path_buff,PATH_MAX)){
      Log::Error<<"Get current absolute path."<<std::endl;
      boost::throw_exception(std::runtime_error("Error is in getting current absolute path!! Check the file: CIOMonitor::CIOMonitor."));//
    }
    m_root_WorkingPath = std::string(abs_path_buff);
}

CIOMonitor::~CIOMonitor()
{
    //nothing
}

void CIOMonitor::createWorkingPath()
{

}
void CIOMonitor::init()
{
    /** \brief  create working environment
    *
    * \path1: scratch
    * \path2: working     create  pop_num  folder
    * \path3: back_up     initial files
    *
    */
}
void CIOMonitor::absoluteFilePath(size_t generation,size_t population,std::string& path_str)
{
    std::stringstream str;

}
std::string& CIOMonitor::currentPath()
{
   return this->m_root_WorkingPath;
}
bool CIOMonitor::checkExeNecessaryFiles(std::vector<std::string>& files, size_t checkmode)
{
	// check if dir_name is a valid dir
	struct stat s;
	lstat(this->m_root_WorkingPath.c_str(), &s );
	if(!S_ISDIR( s.st_mode )){
        Log::Error<<m_root_WorkingPath << "is not valid path."<<std::endl;
        boost::throw_exception(std::runtime_error("Error is in working path!! Check the file: CIOMonitor::checkExeNecessaryFiles"));//
	}

	struct dirent * filename;    // return value for readdir()
 	DIR * dir;                   // return value for opendir()
	dir = opendir(m_root_WorkingPath.c_str());
	if( NULL == dir ){
		Log::Error<<"Can not open dir "<<m_root_WorkingPath<<std::endl;
		boost::throw_exception(std::runtime_error("Error is in working path!! Check the file: CIOMonitor::checkExeNecessaryFiles"));//
	}

	/* read all the files in the dir ~ */
	while( ( filename = readdir(dir) ) != NULL )
	{
		// get rid of "." and ".."
		if( strcmp( filename->d_name , "." ) == 0 ||
			strcmp( filename->d_name , "..") == 0)
			continue;
        lstat(filename->d_name,&s);
        if(S_ISREG(s.st_mode)){

        }
	}
}

}//namespace
