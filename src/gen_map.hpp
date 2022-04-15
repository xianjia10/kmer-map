#include<string>

//read file and init reference file index and reads file index
int read_file(std::string paffiles,std::string reffile,std::string readsfile);
//genrate kmer pos map by paffile
int generate_posmap(std::string b[],FILE *fp,FILE *fw);