#include<bits/stdc++.h>
using namespace std;

typedef pair<uint32_t,uint64_t> KMER;
typedef vector<KMER> KMER_LIST;
typedef unordered_map<string,pair<uint32_t,uint32_t>> lib;
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
    uint32_t fposstart,fposend;
    int flag=0;
    ifstream  posfile(path);
    while(getline(posfile,buf))
    {
        if(buf[0]=='@'){
        if(flag==0)
        {
            name=buf.substr(1,buf.size());
            fposstart=posfile.tellg();
            flag=1;
         } 
        else
        {
            library[name]={fposstart,fposend};
            name=buf.substr(1,buf.size());
            fposstart=posfile.tellg();
        }
        }
        fposend=posfile.tellg();
    }
    library[name]={fposstart,fposend};
    cout<<"library success"<<endl;
}

void ref_lib_init(string name ,string path)
{
    uint32_t p,fposstart,fposend;
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

unordered_map<uint64_t,uint32_t> searchinref(int startpos,int endpos)
{
    unordered_map<uint64_t,uint32_t> klist;
    auto iters= lower_bound(ref_data.begin(),ref_data.end(),startpos,compare);
    //cout<<iters->first<<endl;
    auto itere = lower_bound(ref_data.begin(),ref_data.end(),endpos,compare);
    while(iters!=itere)
    {
        klist[iters->second]=iters->first;
        iters++;
    }
    return klist;
}


vector<pair<uint64_t,uint32_t>> searchinreads(string name,int startpos,int endpos,FILE *fp)
{
    pair<uint64_t,uint32_t> pos;
    vector<pair<uint64_t,uint32_t>> klist;
    string buf;
    uint32_t fposstart,fposend,p;
    uint64_t k;
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

void searchbypaf(string paffiles,string reffile,string readsfile)
{
    string temp;
    string name=" ";
    string b[100];
    unordered_map<uint64_t,uint32_t> ref_list;
    vector<pair<uint64_t,uint32_t>> reads_list;
    ifstream paffile(paffiles);
    lib_init(reffile,ref_lib);
    lib_init(readsfile,reads_lib);
    FILE *fw=NULL,*fp=NULL;
    string path=paffiles+"regen.pos";
    fw=fopen(path.c_str(),"w");
    if (fw == NULL)
    {
        perror("file fopen error!");
        return;
    }
    fp=fopen(readsfile.c_str(),"r");
    if (fp == NULL)
    {
        perror("file fopen error!");
        return;
    }
    while(getline(paffile,temp))
    {
        //cout<<temp<<endl;
        
        split(temp,b);
        //cout<<"read info: "<<b[5]<<atoi(b[7].c_str())<<atoi(b[8].c_str())<<endl;
        if(name!=b[5])
        {
            ref_data.clear();
            ref_lib_init(b[5],reffile);
            name=b[5];
        }
        ref_list=searchinref(atoi(b[7].c_str()),atoi(b[8].c_str()));
        //cout<<"searchinref "<<ref_list.size()<<endl;
        reads_list=searchinreads (b[0],atoi(b[2].c_str()),atoi(b[3].c_str()),fp);
        //cout<<"searchinreads "<<reads_list.size()<<endl;
        for(int i=0;i<9;i++)
        {
            fprintf(fw, "%s\t",b[i].c_str());
        }
        auto iter_reads=reads_list.begin();
        while(iter_reads!=reads_list.end())
        {
            auto ietr_ref=ref_list.find(iter_reads->first);
            if(ietr_ref!=ref_list.end())
            {
                fprintf(fw, "%d,%d\t", iter_reads->second,ietr_ref->second);
            }
            iter_reads++;
        }
        fprintf(fw, "\n");

    }
    fclose(fw);
    fclose(fp);
}

int main(int argc,char **argv)
{
    if(argc<4)
    {
        cout<<"too few argv:1.paf 2.referencepos 3.readspos"<<endl;
        exit(0);
    }
    string paffiles,reffile,readsfile;
    paffiles=argv[1];
    reffile=argv[2];
    readsfile=argv[3];
    searchbypaf(paffiles,reffile,readsfile);
    return 0;
}
