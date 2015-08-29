#include <cstdio>
#include <cstdlib>
#include <map>
#include <cstring>
#include <assert.h>
struct cmp_str { 
  bool operator()(const char *a,const  char *b) { 
    return std::strcmp(a, b) < 0; 
  } 
}; 


typedef std::map<char*,char *,cmp_str> asso;

#define LENS 10000
char buf[LENS];
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


int main(int argc,char**argv){
  
  FILE *fp = fopen(argv[1],"rb");
  asso en;
  while(fgets(buf,LENS,fp)){
    char *tok = strtok(buf,"\n\r\t ");
    char *rs=strtok(NULL,"\n\r\t ");
    char *pos1=strtok(NULL,"\n\r\t ");
    char *pos=strtok(NULL,"\n\r\t ");
    char *key=(char*)calloc(LENS,sizeof(char));
    snprintf(key,LENS,"%s_%s",tok,pos);
    assert(en.find(key)==en.end());
    en[key]=strdup(rs);
  }
  fprintf(stderr,"nitems in bim1:%lu\n",en.size());
  
  fclose(fp);
  fp = fopen(argv[2],"rb");
  int rsid=0;
  while(fgets(buf,LENS,fp)){
    char *tok = strtok(buf,"\n\r\t ");
    char *rs=strtok(NULL,"\n\r\t ");
    char *pos1=strtok(NULL,"\n\r\t ");
    char *pos=strtok(NULL,"\n\r\t ");
    char *al1=strtok(NULL,"\n\r\t ");
    char *al2=strtok(NULL,"\n\r\t ");
    char *key=(char*)calloc(LENS,sizeof(char));
    snprintf(key,LENS,"%s_%s",tok,pos);
    asso::iterator it = en.find(key);
    char b1=tob(al1);
    char b2=tob(al2);
    if(it==en.end()){
      fprintf(stdout,"%s\trs%d\t%s\t%s\t%c\t%c\n",tok,rsid++,pos1,pos,b1,b2);
    }else{
      fprintf(stdout,"%s\t%s\t%s\t%s\t%c\t%c\n",tok,it->second,pos1,pos,b1,b2);
    }
  }
  fclose(fp);
  return 0;
}
