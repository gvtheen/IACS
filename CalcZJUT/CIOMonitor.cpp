#include <stdio.h>
#include <limits.h>
#include "unistd.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "CIOMonitor.h"
#include "../Util/log.hpp"
#include "../GACatalyst.h"

using util::Log;

namespace CALCZJUT{

CIOMonitor::CIOMonitor(CParameter* m_parameter)
:m_pParameter(m_parameter)
{
    //ctor
}

CIOMonitor::~CIOMonitor()
{
    //nothing
}

void CIOMonitor::createWorkingPath()
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
    str<<m_currentPath<<"/Gen_"<<generation<<"_Pop_"<<population;
    str>>path_str;
}
std::string& CIOMonitor::currentPath()
{
   char abs_path_buff[PATH_MAX];
   if(getcwd(abs_path_buff,PATH_MAX)==NULL)
   {
      Log::Error<<"Error in getting current path. CIOMonitor::currentPath!"<<std::endl;
      boost::throw_exception(std::runtime_error("Error in getting current path. CIOMonitor::currentPath!"));
   }
   m_currentPath = abs_path_buff;
   return m_currentPath;
}


}//namespace
