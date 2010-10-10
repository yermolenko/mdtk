#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <mdtk/SimLoop.hpp>


int
main(int argc , char *argv[])
{
  /*
  std::stringstream ss_txt(std::stringstream::in | std::stringstream::out);
  yaatk::unzip_stringstream("in.txt.gz",ss_txt);
  TRACE("Unzipped Ok!");
  yaatk::zip_stringstream("out.txt.gz",ss_txt);
  TRACE("Zipped Ok!");

  std::stringstream ss1_bin(std::stringstream::binary | std::stringstream::in | std::stringstream::out);
  yaatk::unzip_stringstream("in.bin.gz",ss1_bin);
  TRACE("Unzipped Ok!");
  yaatk::zip_stringstream("out.bin.gz",ss1_bin);
  TRACE("Zipped Ok!");
  */


  {
  yaatk::Stream::zipCmdGlobal="xz";

  yaatk::Stream os("dummy.txt",true,false);
  os << "Test!!!";
  os.close();

  yaatk::Stream is("dummy.txt",false,false);
  std::string s;
  is >> s;
  std::cout << s << "\n";
  is.close();
  }
  
  {
  yaatk::Stream::zipCmdGlobal="gzip";

  yaatk::Stream s("in.txt",false,false);
  s.setOutputMode(true);
  s.setFileName("out.txt");
  s.close();
  }

  {
  yaatk::Stream::zipCmdGlobal="gzip";

  yaatk::Stream s("in.bin",false,true);
  s.setOutputMode(true);
  s.setFileName("out.bin");
  s.close();
  }
  
  {
  yaatk::Stream::zipCmdGlobal="gzip";
  yaatk::binary_ifstream s("in.bin");
  s.setOutputMode(true);
  s.setFileName("out.bis.bin");
  s.zipCmd="bzip2";
  s.close();
  }

  {
  yaatk::Stream::zipCmdGlobal="gzip_internal";
  yaatk::binary_ifstream s("in.bin");
  s.setOutputMode(true);
  s.setFileName("out.bis_internal.bin");
  s.close();
  }

  {
  yaatk::Stream::zipCmdGlobal="gzip_internal";
  yaatk::binary_ifstream s("in.bin.gz");
  s.setOutputMode(true);
  s.setFileName("out.bin.guessed");
  s.close();
  }

  {
  yaatk::Stream::zipCmdGlobal="nozip";
  yaatk::binary_ifstream s("test1.cxx");
  s.setOutputMode(true);
  s.setFileName("test1.cxx.yaatk");
  s.close();
  }

  {
  yaatk::Stream::zipCmdGlobal="gzip_dontexist";
  yaatk::binary_ifstream s("in.bin.de");
  s.setOutputMode(true);
  s.setFileName("out.bis_internal.bin.de");
  s.close();
  }
  
  return 0;
}
