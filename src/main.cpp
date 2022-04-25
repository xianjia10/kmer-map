#include "clipp.h"
#include <iostream>
#include "kmerfind.hpp"
#include "generatemap.hpp"
#include <thread>
#include<unistd.h>


using namespace clipp;using std::cout;using std::string;using std::thread;
/*
	description:



*/

string program_name="kmermap";
string version="V0.0.1";


int main(int argc, char* argv[]) {
	bool help=false;
	char mode = 'N';
    string kmer_file = "";
    string ref_file = "";
	string query_file = "";
	string output_path = "./";
	string align = "";
	string rpos="",qpos="";
	bool align_fmt = 0;
    int t1=3,t2=8;
	auto build = "build:" %(
		option("-q", "--query").doc("Query/Reads fasta file") & value("query",query_file),
		option("-r", "--reference").doc("Reference fasta file") & value("reffile", ref_file),
		value("kmerfile", kmer_file).doc("Kmer file with creating by jellyfish")
	);
	auto gen_map = "genmap:" %(
		value("r_pos",rpos).doc("Ref pos file built with build"),
		value("q_pos",qpos).doc("Query pos file built with build"),
		value("paf", align).doc("Alignment file")
	);


    auto cli = (
		option("-t", "--threads").doc("Number of threads") & value("threads", t2),
		option("-o", "--out-put").doc("Diractory of output") & value("out_path", output_path),
		(command("build").set(mode,'b'),build) | (command("genmap").set(mode,'g'),gen_map) | command("-h","--help").set(help,true).doc("Show this page") | command("-V","--version")([](){ cout<<program_name<<"__"<<version<<"__"<<std::endl;}).doc("Version")
    );

	auto fmt = doc_formatting{}
		.first_column(8)                           //left border column for text body
		.doc_column(30)                            //column where parameter docstring starts
		.last_column(100)                          //right border column for text body
		.indent_size(4)                            //indent of documentation lines for children of a documented group
		.line_spacing(0)                           //number of empty lines after single documentation lines
		.paragraph_spacing(1)                      //number of empty lines before and after paragraphs
		.flag_separator(", ")                      //between flags of the same parameter
		.param_separator(" ")                      //between parameters 
		.group_separator(" ")                      //between groups (in usage)
		.alternative_param_separator("/")          //between alternative flags 
		.alternative_group_separator(" | ")        //between alternative groups 
		.surround_group("(", ")")                  //surround groups with these 
		.surround_alternatives("(", ")")           //surround group of alternatives with these
		.surround_alternative_flags("", "")        //surround alternative flags with these
		.surround_joinable("(", ")")               //surround group of joinable flags with these
		.surround_optional("[", "]")               //surround optional parameters with these
		.surround_repeat("", "...")                //surround repeatable parameters with these
		.empty_label("")                           //used if parameter has no flags and no label
		.max_flags_per_param_in_usage(1)           //max. # of flags per parameter in usage
		.max_flags_per_param_in_doc(32)            //max. # of flags per parameter in detailed documentation
		.split_alternatives(true)                  //split usage into several lines for large alternatives
		.alternatives_min_split_size(3)            //min. # of parameters for separate usage line
		.merge_alternative_flags_with_common_prefix(false)  //-ab(cdxy|xy) instead of -abcdxy|-abxy
		.ignore_newline_chars(false)               //ignore '\n' in docstrings
		;


    if(!parse(argc, const_cast<char **>(argv), cli)) {
		cout << "Usage:\n" << usage_lines(cli, "kmermap", fmt)
     << "\nOptions:\n" << documentation(cli, fmt) << "\nERROR: Required parameter missing\n";
		// throw "Division by zero condition!";
		exit(0);
	};
	if (output_path!="./"){
		string cmd = "mkdir -p"+output_path;
		system(cmd.c_str());
	}
	if (help){
		cout << "Usage:\n" << usage_lines(cli, "kmermap", fmt)
     << "\nOptions:\n" << documentation(cli, fmt) << '\n';
		// throw "Division by zero condition!";
		return 0;
	}

	
	/* 
	 *	Build the kmer pos for certain fasta file
	*/
	if (mode=='b'){
		cout<<"Building kmer pos for input fa\nReading the kmer file----------------------------";
		if (fopen(kmer_file.c_str(),"r")){
			read_kmer(kmer_file.c_str(),t1);
			cout<<"Saved the kmer in hash table\n";
		}
		else{
			perror("kmer_file not exist!");
			exit(0);
		}
		if (fopen(ref_file.c_str(),"r"))
			build_pos(ref_file.c_str(), output_path, 1, t1, t2);
		else{
			perror("Reference file not exist!");
			exit(0);
		}
		if (fopen(query_file.c_str(),"r"))
			build_pos(query_file.c_str(),output_path, 0, t1, t2);
		else{
			perror("Query file not exist!");
			exit(0);
		}
		cout<<"Building pos finished!\n";
	}


	if (mode=='g'){
		if (align != ""){
			// align_fmt will be transmited later.
			string newalign=align.substr(0,align.size()-3)+"srt.paf";
			if(access(newalign.c_str(),F_OK)==-1){
				string cmd="sort -k6,6 -k1,1 "+align+" > "+newalign; 
				system(cmd.c_str());
				cout<<"align file trans"<<std::endl;
			}
			read_file(newalign,rpos,qpos,output_path);
		}
	}
	return 0;
}