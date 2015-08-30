#include <cstdio>
#include <cstdlib>
#include <map>
#include <cstring>
#include <assert.h>
#include "../reffinder/refFinder.h"



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


typedef struct{
  char *rs;
  char a1;
  char a2;
}dat;



typedef std::map<char*,dat,cmp_str> asso;
#define LENS 2048

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



int main(int argc,char**argv){
  if(argc==1){
    fprintf(stderr,"possible options\n\t 1) checkstrand hg19.fa.gz file.bim\n\t\t This will give you the entries that doesnt match \n");
    return 0;
  }
  argc--;++argv;  
  if(!strcasecmp(*argv,"checkstrand"))
    return checkstrand(--argc,++argv);

  FILE *fp = fopen("filtered.bim","rb");//<-assumed hg19 strand and super
  asso en;
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
