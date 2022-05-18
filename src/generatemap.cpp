#include<bits/stdc++.h>
#include"generatemap.hpp"
using namespace std;

typedef pair<uint32_t,uint64_t> KMER;
typedef vector<KMER> KMER_LIST;
typedef unordered_map<string,pair<uint64_t,uint64_t>> lib;
lib ref_lib,reads_lib; 
KMER_LIST ref_data;

void split(string q,string item[])
{
    stringstream ss;
    int i=0;
    ss.clear();
    ss.str(q); 
    while(!ss.fail())
    {
        ss>>item[i];
        i++;
    }
}

//initize hash
void lib_init(string path,lib &library)
{
    string buf,name;
    uint64_t fposstart,fposend;
    int flag=0;
    ifstream  posfile(path);
    while(getline(posfile,buf))
    {
        if(buf[0]=='@'){
            if(flag==0)
            {
                name=buf.substr(1,buf.size()-1);
                fposstart=posfile.tellg();
                flag=1;
            } 
            else
            {
                library[name]={fposstart,fposend};
                name=buf.substr(1,buf.size()-1);
                fposstart=posfile.tellg();
            }
            }
            fposend=posfile.tellg();
    }
    library[name]={fposstart,fposend};
    cout<<"library success"<<endl;
}

//read reference
void ref_lib_init(string name ,string path)
{
    uint32_t p;
    uint64_t fposstart,fposend;
    fposstart=ref_lib[name].first;
    fposend=ref_lib[name].second;
    uint64_t k;
    string buf;
    FILE *fp;
    fp=fopen(path.c_str(),"r");
    fseek(fp,fposstart,SEEK_SET);
    while(ftell(fp)<fposend)
    {
        fscanf(fp,"%8x %11lx",&p,&k);
        ref_data.push_back(make_pair(p,k));
    }
    fclose(fp);
}

inline bool compare(KMER a,uint32_t b)
{
    return a.first<b; 
}

//search kmer in reference
unordered_map<uint64_t,uint32_t> searchinref(int startpos,int endpos)
{
    unordered_map<uint64_t,uint32_t> klist;
    auto iters= lower_bound(ref_data.begin(),ref_data.end(),startpos,compare);
    auto itere = lower_bound(ref_data.begin(),ref_data.end(),endpos,compare);
    while(iters!=itere)
    {
        klist[iters->second]=iters->first;
        iters++;
    }
    return klist;
}

//search kemr in query 
vector<pair<uint64_t,uint32_t>> searchinreads(string name,int startpos,int endpos,FILE *fp)
{
    pair<uint64_t,uint32_t> pos;
    vector<pair<uint64_t,uint32_t>> klist;
    string buf;
    uint32_t p;
    uint64_t k,fposstart,fposend;
    int num;
    fposstart=reads_lib[name].first;
    fposend=reads_lib[name].second;
    num=(fposend - fposstart)/21;
    int lo = 0, hi = num;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        fseek(fp,fposstart+mid*21,SEEK_SET);
        fscanf(fp,"%8x%11lx",&p,&k);
        if (p < startpos) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    fseek(fp,fposstart+lo*21,SEEK_SET);
    while(ftell(fp)<fposend)
    {
        fscanf(fp,"%8x%11lx",&p,&k);
        if(p>endpos)    break;
        klist.push_back(make_pair(k,p));
    }
    return klist;
}

//genrate kmer pos map
int generate_posmap(string b[],FILE *fp,FILE *fw)
{
    unordered_map<uint64_t,uint32_t> ref_list;
    vector<pair<uint64_t,uint32_t>> reads_list;
    unordered_map<uint64_t,uint32_t> ref_list_over;
    int kmapq;
    int common_kmer=0;
    int count=0;
    ref_list=searchinref(atoi(b[7].c_str()),atoi(b[8].c_str()));
    reads_list=searchinreads (b[0],atoi(b[2].c_str()),atoi(b[3].c_str()),fp);
    int fenmu=min(ref_list.size(),reads_list.size());
    for(int i=0;i<9;i++)
    {
        fprintf(fw, "%s\t",b[i].c_str());
    }
    auto iter_reads=reads_list.begin();
    while(iter_reads!=reads_list.end())
    {
        auto ietr_ref=ref_list.find(iter_reads->first);
        auto iter_ref_over=ref_list_over.find(iter_reads->first);
        if(ietr_ref!=ref_list.end())
        {
            fprintf(fw, "%d,%d\t", iter_reads->second,ietr_ref->second);
            common_kmer++;
            ref_list_over[ietr_ref->first]=ietr_ref->second;
            ref_list.erase(ietr_ref->first);
        }
        else if(iter_ref_over!= ref_list_over.end()){
            fprintf(fw, "%d,%d\t", iter_reads->second,iter_ref_over->second);
        }
        else{
            fprintf(fw, "%d,%d\t", iter_reads->second,0);
        }
        count++;
        if(count==10000)
        {
            fprintf(fw,"-1\n");
            for(int i=0;i<9;i++)
            {
                fprintf(fw, "%s\t",b[i].c_str());
            }
            count=0;
        }
        iter_reads++;
    }
    for(auto it : ref_list){
	    fprintf(fw, "%d,%d\t",0,it.second);
        count++;
        if(count==10000)
        {
            fprintf(fw,"-1\n");
            for(int i=0;i<9;i++)
            {
                fprintf(fw, "%s\t",b[i].c_str());
            }
            count=0;
        }
    }
    if (fenmu==0)
    {
        kmapq=0;
    }
    else{
        float v=(float)(common_kmer)/fenmu;
        kmapq=(int)(100*v+0.5);
    }
    fprintf(fw, "%d\n",kmapq);
    return 1;
}

//read paf file and init data's index
int read_file(string paffiles,string reffile,string readsfile,string out_path)
{
    string temp;
    string name="*";
    string b[100];
    uint64_t fpos;
    ifstream paffile(paffiles);
    ref_lib.clear();
    reads_lib.clear();
    lib_init(reffile,ref_lib);
    lib_init(readsfile,reads_lib);
    FILE *fw=NULL,*fp=NULL,*ri=NULL,*qi=NULL;
    string path=out_path+"/"+paffiles.substr(0,paffiles.size()-3)+"kmermap";
    fw=fopen(path.c_str(),"w");
    if (fw == NULL){
        perror("file fopen error!");
        exit(0);
    }
    fp=fopen(readsfile.c_str(),"r");
    if (fp == NULL){
        perror("file fopen error!");
        exit(0);
    }
    string pathr=path+".reference.index";
    string pathq=path+".query.index";
    ri=fopen(pathr.c_str(),"w");
    qi=fopen(pathq.c_str(),"w");
    if(ri==NULL||qi==NULL)
    {
        perror("file fopen error!");
        exit(0);
    }
    while(getline(paffile,temp)){ 
        split(temp,b);
        fpos=ftell(fw);
        if(name!=b[5]){
            ref_data.clear();
            ref_lib_init(b[5],reffile);
            name=b[5];
            fprintf(ri,"%s\t%ld\n",name.c_str(),fpos);
        }
        fprintf(qi,"%s\t%ld\n",b[0].c_str(),fpos);
        generate_posmap(b,fp,fw);
    }
    fclose(fw);
    fclose(fp);
    fclose(ri);
    fclose(qi);
    return 1;
}
