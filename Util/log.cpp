/******************************************************************************
**
** Copyright (C) 2019-2031 Dr.Gui-lin Zhuang <glzhuang@zjut.edu.cn>
** All rights reserved.
**
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/
/**
 * @file log.cpp
 * @author Matthew Amidon
 *
 * Implementation of the Log class.
 *
 * mlpack is free software; you may redistribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#include "log.hpp"


using namespace util;
//std::ostream& Log::cout ;

std::ofstream Log::fcout("IACS.log",std::ios::out);

PrefixedOutStream Log::Warn(Log::fcout,"[Warn]: ");
//PrefixedOutStream Log::Warn(std::cout,"[Warn]: ");

PrefixedOutStream Log::Info(Log::fcout,"[Info]: ");
//PrefixedOutStream Log::Info(std::cout,"[Info]: ");

#ifdef DEBUG
  PrefixedOutStream Log::Debug(Log::fcout,"[Debug]: ");
#else
  NullOutStream Log::Debug;
#endif

PrefixedOutStream Log::Fatal(Log::fcout,"[Fatal]: ");

PrefixedOutStream Log::Error(Log::fcout,"[Error]: ");

PrefixedOutStream Log::Output(Log::fcout,"");
// Only do anything for Assert() if in debugging mode.


#ifdef DEBUG
void Log::Assert(bool condition, const std::string& message)
{
  if (!condition)
  {
#ifdef HAS_BFD_DL
    Backtrace bt;

    Log::Debug << bt.ToString();
#endif
    Log::Debug << message << std::endl;

    throw std::runtime_error("Log::Assert() failed: " + message);
  }
}
#else
void Log::Assert(bool /* condition */, const std::string& /* message */)
{ }


#endif
