#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <mdtk/SimLoop.hpp>
#include <sstream>


int
main(int argc , char *argv[])
{
  //  std::stringstream ss;
  std::stringstream ss_txt;
  yaatk::unzip_stringstream("in.txt.gz",ss_txt);
  TRACE("Unzipped Ok!");
  yaatk::zip_stringstream("out.txt.gz",ss_txt);
  TRACE("Zipped Ok!");

  std::stringstream ss1_bin;//(std::stringstream::binary);
  yaatk::unzip_stringstream("in.bin.gz",ss1_bin);
  TRACE("Unzipped Ok!");
  yaatk::zip_stringstream("out.bin.gz",ss1_bin);
  TRACE("Zipped Ok!");

  return 0;
}
