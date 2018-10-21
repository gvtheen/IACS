#include "CProgameDoc.h"
#include "log.hpp"

namespace util{

CProgameDoc::CProgameDoc()
{
    this->versionText = " 0.1 ";
    this->authorText = "Author: Gui-lin Zhuang ( glzhuang@zjut.edu.cn )";
    this->addressText= "Institute: Zhejiang University of Technology";
}
void CProgameDoc::output()
{
    Log::Output<<"******************************************************************************"<<std::endl;
    Log::Output<<"******************************************************************************"<<std::endl;
    Log::Output<<"**   Copyright (C) 2019-2031 Dr.Gui-lin Zhuang <glzhuang@zjut.edu.cn>       **"<<std::endl;
    Log::Output<<"**                           Dr.Jian-guo Wang  <jgw@zjut.edu.cn>            **"<<std::endl;
    Log::Output<<"**     ^^^^^^          ^^^^             ^^^^^^           ^^^^^^^            **"<<std::endl;
    Log::Output<<"**     ^^^^^^         ^^^^^^          ^^^^^^^^^         ^^^^^^^^^           **"<<std::endl;
    Log::Output<<"**       ^^           ^^  ^^          ^^                ^^     ^^           **"<<std::endl;
    Log::Output<<"**       ^^           ^^  ^^          ^^                 ^^                 **"<<std::endl;
    Log::Output<<"**       ^^           ^^  ^^          ^^                  ^^                **"<<std::endl;
    Log::Output<<"**       ^^           ^^^^^^          ^^                   ^^               **"<<std::endl;
    Log::Output<<"**       ^^           ^^^^^^          ^^                    ^^              **"<<std::endl;
    Log::Output<<"**       ^^           ^^  ^^          ^^                ^^   ^^             **"<<std::endl;
    Log::Output<<"**     ^^^^^^         ^^  ^^          ^^^^^^^^^         ^^^^^^^^            **"<<std::endl;
    Log::Output<<"**     ^^^^^^         ^^  ^^            ^^^^^^           ^^^^^^             **"<<std::endl;
    Log::Output<<"**                                                                          **"<<std::endl;
    Log::Output<<"******************************************************************************"<<std::endl;
    Log::Output<<"******************************************************************************"<<std::endl;
    Log::Output<<"Current Version: IACS" <<this->versionText<<std::endl;
    Log::Output<<this->authorText <<";" <<this->addressText<<std::endl;
}
CProgameDoc::~CProgameDoc()
{
    //dtor
}

}
