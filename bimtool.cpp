#include <cstdio>
#include <cstdlib>
#include <map>
#include <cstring>
#include <assert.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>
#include <fcntl.h>

#include "../reffinder/refFinder.h"

#define LENS 2048

int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}

int copy_file (char *in,char *out){
  fprintf(stderr,"\t-> Will copy: %s to %s\n",in,out);
  //  int BUFSIZ = 8192;
  char buf[BUFSIZ];
  size_t size;
  
  int source = open(in, O_RDONLY, 0);
  int dest = open(out, O_WRONLY | O_CREAT /*| O_TRUNC/**/, 0644);
  
  while ((size = read(source, buf, BUFSIZ)) > 0) {
    write(dest, buf, size);
  }
  
  close(source);
  close(dest);
  

  return 0;
}

char *makename(char *name,char *extension){
  char *newname=new char[LENS];
  newname=strcat(newname,name);
  newname=strcat(newname,extension);
  while(fexists(newname))
    newname=strcat(newname,extension);
  
  fprintf(stderr,"\t-> Output filename will be: %s\n",newname);
  fprintf(stderr,"%s\n",name);
  return newname;
}



struct cmp_str { 
  bool operator()(const char *a,const  char *b) { 
    return std::strcmp(a, b) < 0; 
  }
}; 

int refToInt[256] = {
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};

char intToRef[5] = {'A','C','G','T','N'};

struct info{
  char *rs;
  char *chr;
  char *bp;
  int strand;
  char al1;
  char al2;
};

void print(FILE *fp,const info &i){
  fprintf(fp,"rs:%s chr:%s bp:%s strand:%d al1:%c al2:%c\n",i.rs,i.chr,i.bp,i.strand,i.al1,i.al2);

}

struct dat{
  char *rs;
  char a1;
  char a2;
  int strand;//-1 negative; 1 positive; 0 undefined
};

bool operator==(const info &lhs,const info &rhs){

  if(strcmp(lhs.rs,rhs.rs)){
    //  fprintf(stderr,"lhs:%s rhs:%s\n",lhs.rs,rhs.rs);
    return false;
  }if(strcmp(lhs.chr,rhs.chr))
    return false;
  if(lhs.strand!=rhs.strand)
    return false;
  if(lhs.al1!=rhs.al1)
    return false;
  if(lhs.al2!=rhs.al2)
    return false;
  return true;
}
bool operator!=(const info &lhs,const info &rhs){
  return !(lhs==rhs);
}




char *buf=new char[LENS];
char *buf2=new char[LENS];
#include <ctype.h>
char tob(char *c){
  int i=atoi(c);
  if(i==1)
    return 'A';
  if(i==2)
    return 'C';
  if(i==3)
    return 'G';
  if(i==4)
    return 'T';
  return '0';
}

int flip(char c){
  if(c=='A')
    return 'T';
  if(c=='C')
    return 'G';
  if(c=='G')
    return 'C';
  if(c=='T')
    return 'A';
  return c;
}




int checkstrand(int argc,char**argv){
  fprintf(stderr,"argc:%d argv:%s\n",argc,*argv);
  if(argc!=2){
    fprintf(stderr," checkstrand hg19.fa.gz file.bim\n");
    return 0;
  }else{
    fprintf(stderr,"\t assuming fasta=%s bimfile=%s\n",argv[0],argv[1]);

  }
  fprintf(stderr,"\t-> Program will only use sites where both alleles are observed in data\n");
  perFasta *ref = init(argv[0]);
  
  size_t mismatch[5][5][5];
  for(int i=0;i<5;i++)
    for(int ii=0;ii<5;ii++)
      for(int iii=0;iii<5;iii++)
	mismatch[i][ii][iii]=0;

  FILE *fp = fopen(argv[1],"rb");//<-assumed hg19 strand and super
  while(fgets(buf,LENS,fp)){
    buf2=strcpy(buf2,buf);
    char *tok = strtok(buf,"\n\r\t ");
    char *rs=strtok(NULL,"\n\r\t ");
    char *pos1=strtok(NULL,"\n\r\t ");
    char *pos=strtok(NULL,"\n\r\t ");
    char al1=strtok(NULL,"\n\r\t ")[0];
    char al2=strtok(NULL,"\n\r\t ")[0];
    if(al1=='0'||al2=='0')
      continue;
    char refchar = toupper(getchar(tok,atoi(pos)-1,ref));      
    int a,b,c;
    a=refToInt[refchar];
    b=refToInt[al1];
    c=refToInt[al2];
    mismatch[a][b][c]++;
    if(a!=b&&a!=c)
      fprintf(stdout,"%s",buf2);
  }
  fprintf(stderr,"hg19.fa allele1 allele2\n");
  size_t dd[2] = {0,0};
  for(int i=0;i<4;i++)
    for(int ii=0;ii<4;ii++)
      for(int iii=0;iii<4;iii++){
	if(mismatch[i][ii][iii])
	  fprintf(stderr,"%c %c %c: %lu\n",intToRef[i],intToRef[ii],intToRef[iii],mismatch[i][ii][iii]);
	if(i!=ii&&i!=iii)
	  dd[0] += mismatch[i][ii][iii];
	else
	  dd[1] += mismatch[i][ii][iii];
      }
  fprintf(stderr,"WRONG:%lu RIGHT:%lu\n",dd[0],dd[1]);
}


int noannoying(int argc,char**argv){
  fprintf(stderr,"argc:%d argv:%s\n",argc,*argv);
  if(argc!=1){
    fprintf(stderr," checkstrandfile.bim\n");
    return 0;
  }else{
    fprintf(stderr,"\t assuming bimfile=%s\n",argv[0],argv[1]);
  }
  fprintf(stderr,"\t-> Program will only only print out the variable sites with A<->C and T<->G mutations removed\n");

  FILE *fp = fopen(argv[0],"rb");
  size_t nbad=0;
  while(fgets(buf,LENS,fp)){
    buf2=strcpy(buf2,buf);
    char *tok = strtok(buf,"\n\r\t ");
    char *rs=strtok(NULL,"\n\r\t ");
    char *pos1=strtok(NULL,"\n\r\t ");
    char *pos=strtok(NULL,"\n\r\t ");
    char al1=strtok(NULL,"\n\r\t ")[0];
    char al2=strtok(NULL,"\n\r\t ")[0];
    if(al1=='0'||al2=='0')
      continue;
    int a,b;
    a=refToInt[al1];
    b=refToInt[al2];
    assert(a!=b);
    if(a==1&&b==2||a==2&&b==1){
      nbad++;
      continue;
    }
    if(a==0&&b==3||a==3&&b==0){
      nbad++;
      continue;
    }
    fprintf(stdout,"%s",buf2);
  }
  fprintf(stderr,"\t-> Removed: %lu of the bad annoying snps\n",nbad);
  fclose(fp);
}


int flipstrand(int argc,char**argv){
  fprintf(stderr,"argc:%d argv:%s\n",argc,*argv);
  if(argc!=2){
    fprintf(stderr," flipstrand info file.bim\n");
    return 0;
  }else{
    fprintf(stderr,"\t assuming info=%s bimfile=%s\n",argv[0],argv[1]);

  }
    
  char *newname = makename(argv[1],".original");
  copy_file(argv[1],newname);

  typedef std::map<char*,info,cmp_str> assoRs;
  typedef std::map<char*,info,cmp_str> assoChrPos;
  assoRs as;
  assoChrPos acp;
  std::vector<char *> rsdup;
  std::vector<char *> cpdup;
  FILE *fp = fopen(argv[0],"rb");//<-assumed hg19 strand and super
  while(fgets(buf,LENS,fp)){
    char *rs=strdup(strtok(buf,"\n\r\t "));
    char *chr=strtok(NULL,"\n\r\t ");
    char *pos=strtok(NULL,"\n\r\t ");
    char *strand=strtok(NULL,"\n\r\t ");
    char al1=strtok(NULL,"\n\r\t ")[0];
    char al2=strtok(NULL,"\n\r\t ")[0];
    char *key1 = new char[LENS];
    snprintf(key1,LENS,"%s_%s",chr,pos);
    info in;
    in.rs=strdup(rs);
    in.chr=strdup(chr);
    in.bp=strdup(pos);
    if(strlen(strand)==3)
      in.strand=0;
    else if(strand[0]=='-')
      in.strand=-1;
    else if(strand[0]=='+')
      in.strand=1;
    else{
      fprintf(stderr,"YOGOGOGOGOGOG\n");
    }
    in.al1=al1;
    in.al2=al2;
    //valdiate that duplicate rsnumbers are identical and that duplicate chr_pos are identical
    assoRs::iterator it = as.find(rs);
    int skip=0;
    if(it!=as.end()){
      if(in!=it->second){
	fprintf(stderr,"\t-> Found duplicate rsnumber: %s \n",rs);
	print(stderr,in);
	print(stderr,it->second);
	rsdup.push_back(rs);
      }
    }
    assoChrPos::iterator it2 = acp.find(key1);
    if(it2!=acp.end()){
      if(in!=it2->second){
	fprintf(stderr,"Found duplicate position: %s \n",key1);
	print(stderr,in);
	print(stderr,it2->second);
	cpdup.push_back(key1);
      }
    }
    if(strcmp(rs,"---"))
      as[rs] = in;
    if(strcmp(key1,"---_---"))
      acp[key1] = in;
  }
  fprintf(stderr,"Will remove duplicate rsnumber and chr_pos\n");
  for(unsigned i=0;i<rsdup.size();i++)
    as.erase(rsdup[i]);
  for(unsigned i=0;i<cpdup.size();i++)
    acp.erase(cpdup[i]);
  fclose(fp);
  fprintf(stderr,"\t-> Unique positions in annotaionfile using rsnumbers:%lu using chr_pos:%lu\n",as.size(),acp.size());
  
  fp=fopen(newname,"rb");
  while(fgets(buf,LENS,fp)){
    char *tok = strtok(buf,"\n\r\t ");
    char *rs=strtok(NULL,"\n\r\t ");
    char *pos1=strtok(NULL,"\n\r\t ");
    char *pos=strtok(NULL,"\n\r\t ");
    char al1=strtok(NULL,"\n\r\t ")[0];
    char al2=strtok(NULL,"\n\r\t ")[0];
    char *key=(char*)calloc(LENS,sizeof(char));
    snprintf(key,LENS,"%s_%s",tok,pos);
    assoRs::iterator it = as.find(rs);
    assoChrPos::iterator it2 = acp.find(key);
    info inf;
    if(!strcmp(rs,"---")||!strcmp(key,"---_--=")||!strcmp(key,"0_0"))
      fprintf(stdout,"0\t---\t%s\t0\t%c\t%c\n",pos1,al1,al2);
    else{
      if(it==as.end() && it2==acp.end()){
	fprintf(stdout,"0\t---\t%s\t0\t%c\t%c\n",pos1,al1,al2);
	continue;
      }else if (it!=as.end() && it2!=acp.end()){ 
	if(it->second!=it2->second){
	  fprintf(stderr,"rs number and chr_pos doesnt match :%s vs %s\n",rs,key);
	  fprintf(stdout,"0\t---\t%s\t0\t%c\t%c\n",pos1,al1,al2);
	  continue;
	}else
	  inf=it->second;
      }else if (it!=as.end() && it2==acp.end()) {
	inf=it->second;
      }else if (it==as.end() && it2!=acp.end()) 
	inf=it2->second;
      if(inf.strand==0){
	fprintf(stdout,"0\t---\t%s\t0\t%c\t%c\n",pos1,al1,al2);
	continue;
      }
      fprintf(stdout,"%s\t%s\t%s\t%s\t",inf.chr,inf.rs,pos1,inf.bp);
      if(inf.strand==1)
	fprintf(stdout,"%c\t%c\n",al1,al2);
      else
	fprintf(stdout,"%c\t%c\n",flip(al1),flip(al2));
    }
    
  }

  
  
  return 0;
}



int main(int argc,char**argv){
  if(argc==1){
    fprintf(stderr,"possible options\n\t 1) checkstrand hg19.fa.gz file.bim\n\t\t This will give you the entries that doesnt match \n");
    fprintf(stderr,"\t 2) flipstrand infofile bimfile\n");
    fprintf(stderr,"\t 2) noannoying bimfile\n");
    return 0;
  }
  argc--;++argv;  
  if(!strcasecmp(*argv,"checkstrand"))
    return checkstrand(--argc,++argv);
  if(!strcasecmp(*argv,"flipstrand"))
    return flipstrand(--argc,++argv);
  if(!strcasecmp(*argv,"noannoying"))
    return noannoying(--argc,++argv);
  

  typedef std::map<char*,dat,cmp_str> asso;
  FILE *fp = fopen("filtered.bim","rb");//<-assumed hg19 strand and super
  asso en;
  char buf[LENS];
  while(fgets(buf,LENS,fp)){
    char *tok = strtok(buf,"\n\r\t ");
    char *rs=strtok(NULL,"\n\r\t ");
    char *pos1=strtok(NULL,"\n\r\t ");
    char *pos=strtok(NULL,"\n\r\t ");
    char al1=strtok(NULL,"\n\r\t ")[0];
    char al2=strtok(NULL,"\n\r\t ")[0];
    char *key=(char*)calloc(LENS,sizeof(char));
    snprintf(key,LENS,"%s_%s",tok,pos);
    assert(en.find(key)==en.end());
    dat d;d.rs=strdup(rs);
    d.a1=al1;
    d.a2=al2;
    en[key]=d;
  }
  fprintf(stderr,"nitems in bim1:%lu\n",en.size());
  
  fclose(fp);
  fp = fopen("Beall_DiRienzo_Sherpa/Beall_DiRienzo_Sherpa.bim.original","rb");
  int rsid=0;
  while(fgets(buf,LENS,fp)){
   
    char *tok = strtok(buf,"\n\r\t ");
    char *rs=strtok(NULL,"\n\r\t ");
    char *pos1=strtok(NULL,"\n\r\t ");
    char *pos=strtok(NULL,"\n\r\t ");
    char al1=strtok(NULL,"\n\r\t ")[0];
    char al2=strtok(NULL,"\n\r\t ")[0];
    
    char *key=(char*)calloc(LENS,sizeof(char));
    snprintf(key,LENS,"%s_%s",tok,pos);
    asso::iterator it = en.find(key);
    char b1=toupper(al1);
    char b2=toupper(al2);
    //    fprintf(stderr,"b1:%c b2:%c s1:%c s2:%c\n",b1,b2,al1,al2);
    
    
  
    if(it==en.end()){
      fprintf(stdout,"%s\trrs%d\t%s\t%s\t%c\t%c\n",tok,rsid++,pos1,pos,b1,b2);
    }else{
      fprintf(stdout,"%s\t%s\t%s\t%s\t%c\t%c\n",tok,it->second.rs,pos1,pos,b1,b2);
    }
  }
  fclose(fp);
  return 0;
}
