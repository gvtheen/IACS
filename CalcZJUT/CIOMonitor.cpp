#include <stdio.h>
#include <limits.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include<boost/filesystem.hpp>
#include "CIOMonitor.h"
#include "../Util/log.hpp"
#include "../GACatalyst.h"

using util::Log;

namespace CALCZJUT{

CIOMonitor::CIOMonitor(CParameter* m_parameter)
:m_pParameter(m_parameter)
{
    m_root_WorkingPath = boost::filesystem::current_path();
}

CIOMonitor::~CIOMonitor()
{
    //nothing
}

void CIOMonitor::initWorkEnvironment()
{
    /** \brief  create working environment
    * \path1: scratch
    * \path2: work          create  pop_num  folder
    * \path3: parameter     initial files
    */
    boost::filesystem::path working_path,scratch_path,parameter_path;
      working_path = m_root_WorkingPath/"work";
    if(!boost::filesystem::exists(working_path))
        boost::filesystem::create_directory(working_path);

      scratch_path = m_root_WorkingPath/"scratch";
    if(!boost::filesystem::exists(scratch_path))
        boost::filesystem::create_directory(scratch_path);

    parameter_path = m_root_WorkingPath/"parameter";
    if(!boost::filesystem::exists(parameter_path))
        boost::filesystem::create_directory(parameter_path);

}
void CIOMonitor::setCurrentWorkPathAt(size_t generation,size_t population)
{
    std::stringstream str;
    str<<"Gen_"<<generation<<"Pop_"<<population;
    std::string pathstr;
    str>>pathstr;
    boost::filesystem::path working_path;
    working_path = m_root_WorkingPath/"work";

    working_path=working_path/pathstr.c_str();

    if(! boost::filesystem::exists(working_path))
         boost::filesystem::create_directory(working_path);

    boost::filesystem::current_path(working_path);
}
std::string CIOMonitor::currentInitPath()
{
    return this->m_root_WorkingPath.string();
}
std::string CIOMonitor::currentWorkPath()
{
    return boost::filesystem::current_path().string();
}
void CIOMonitor::checkExeNecessaryFiles(std::vector<std::string>& files)
{
    boost::filesystem::path tempPath;
    for(size_t i=0;i<files.size();i++){
       tempPath=m_root_WorkingPath/files[i].c_str();
       if(!boost::filesystem::is_regular_file(tempPath)||
           boost::filesystem::file_size(tempPath) ==0 ){
           Log::Error<<files[i]<<" parameter file is not exist or empty!"<<std::endl;
           boost::throw_exception(std::runtime_error(files[i] +" parameter file is no exist or empty! CIOMonitor::checkExeNecessaryFiles!\n"));
       }else{
           boost::filesystem::path parameter_path = m_root_WorkingPath/"parameter";
           this->moveFileToPath(tempPath.string(),parameter_path.string());
       }
    }
}
void CIOMonitor::moveFileToPath(const std::string& file, const std::string& dir)
{
    boost::filesystem::path oldFilePath(file);

    boost::filesystem::path newDirPath(dir);

    if(! boost::filesystem::exists(newDirPath))
         boost::filesystem::create_directory(newDirPath);

    boost::filesystem::path newFilePath(dir + "/" + oldFilePath.leaf().string());
    if(!boost::filesystem::is_regular_file(newFilePath) ||
       boost::filesystem::file_size(newFilePath) ==0 )
       boost::filesystem::rename(oldFilePath,newFilePath);
}
void CIOMonitor::copyFileToPath(const std::string& file, const std::string& dir)
{
    boost::filesystem::path oldFilePath(file);

    boost::filesystem::path newDirPath(dir);

    if(! boost::filesystem::exists(newDirPath))
         boost::filesystem::create_directory(newDirPath);

    boost::filesystem::path newFilePath(dir + "/" + oldFilePath.leaf().string());
    if(!boost::filesystem::is_regular_file(newFilePath) ||
       boost::filesystem::file_size(newFilePath) ==0 )
       boost::filesystem::copy_file(oldFilePath,newFilePath);
}


}//namespace
