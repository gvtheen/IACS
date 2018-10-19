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
 * @file prefixedoutstream.cpp
 * @author Ryan Curtin
 * @author Matthew Amidon
 *
 */
//#include "prereqs.hpp"

#include "prefixedoutstream.hpp"

using namespace util;

/**
 * These are all necessary because gcc's template mechanism does not seem smart
 * enough to figure out what I want to pass into operator<< without these.  That
 * may not be the actual case, but it works when these is here.
 */

PrefixedOutStream& PrefixedOutStream::operator<<(bool val)
{
  BaseLogic<bool>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(short val)
{
  BaseLogic<short>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(unsigned short val)
{
  BaseLogic<unsigned short>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(int val)
{
  BaseLogic<int>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(unsigned int val)
{
  BaseLogic<unsigned int>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(long val)
{
  BaseLogic<long>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(unsigned long val)
{
  BaseLogic<unsigned long>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(float val)
{
  BaseLogic<float>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(double val)
{
  BaseLogic<double>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(long double val)
{
  BaseLogic<long double>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(void* val)
{
  BaseLogic<void*>(val);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(const char* str)
{
  BaseLogic<const char*>(str);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(std::string& str)
{
  BaseLogic<std::string>(str);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(std::streambuf* sb)
{
  BaseLogic<std::streambuf*>(sb);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(
    std::ostream& (*pf)(std::ostream&))
{
  BaseLogic<std::ostream& (*)(std::ostream&)>(pf);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(std::ios& (*pf)(std::ios&))
{
  BaseLogic<std::ios& (*)(std::ios&)>(pf);
  return *this;
}

PrefixedOutStream& PrefixedOutStream::operator<<(
    std::ios_base& (*pf) (std::ios_base&))
{
  BaseLogic<std::ios_base& (*)(std::ios_base&)>(pf);
  return *this;
}
