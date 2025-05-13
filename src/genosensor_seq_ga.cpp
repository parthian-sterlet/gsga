#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define NOLIG 30 // max number of oligonucleotides = 2*10
#define NMUT 50 //max number of nucleotide in one mutation block per one mutation
#define PWMLEN 32//max length of PWM
#define LENGTH 350 //dlina polimera = 21*10 
#define NPVT 7000 // 4islo porogov zadannyh p-value
#define POOL 200 // size of population = 4islo osobey populyacii
#define ELIT 25 // size of population = 4islo osobey populyacii - na vyhod elita
#define NMATRICES 600

struct town {
	char mp[NOLIG][NMUT];// eti est' = order
	double fit; //fitness
	double bpv_n;// best score negative
	double bpv_wn;//best score negative weighted
	double bpv_wy;// av best score positive
	int bpos;
	int bmat;
	char ori;
	char seq[PWMLEN];
	void init(char **wt_nt, int *n_mut, int n_olig,int *n_mut_deg);
	void init0(char **wt_nt, int *n_mut, int n_olig);
	void get_copy(town *a, int *n_mut, int n_olig);
	void print_all(int *n_mut, int n_olig, FILE *out);	
} dmut, drec[2], pop[POOL];
void town::init0(char **wt_nt, int *n_mut, int n_olig)
{
	int i,j;
	fit=bpv_wn=bpv_wy=bpv_n=0;
    bmat=0, bpos=0, ori='+';	
	memset(seq,'\0',sizeof(seq));	
	for(j=0;j<n_olig;j++)
	{
		for(i=0;i<n_mut[j];i++)mp[j][i]=wt_nt[j][i];
		mp[j][n_mut[j]] = '\0';
	}	
}
void town::init(char **wt_nt, int *n_mut, int n_olig,int *n_mut_deg)
{
	int i,j,k,s;
	fit=bpv_wn=bpv_wy=bpv_n=0;
    bmat=0, bpos=0, ori='+';	
	memset(seq,'\0',sizeof(seq));	
	for(j=0;j<n_olig;j++)
	{
		for(i=0;i<n_mut[j];i++)mp[j][i]=wt_nt[j][i];
		mp[j][n_mut[j]] = '\0';
	}
	char triplet[4][4]={"cgt","agt","act","acg"};
	char muts[LENGTH];
	for(j=0;j<n_olig;j++)
	{
		int n_mut_deg_cur=0;
		for(k=0;k<n_mut_deg[j];k++)
		{
			int gom, r1;
			do
			{
				gom=1;
				r1=rand()%n_mut[j];			
				for(s=0;s<n_mut_deg_cur;s++)
				{
					if(muts[s]==r1)
					{
						gom=0;
						break;
					}
				}
			}
			while(gom==0);
			muts[n_mut_deg_cur]=r1;
			n_mut_deg_cur++;
			int r2=rand()%3;
			switch (mp[j][r1]) {
				case 'a': {mp[j][r1]=triplet[0][r2];break;}
				case 'c': {mp[j][r1]=triplet[1][r2];break;}
				case 'g': {mp[j][r1]=triplet[2][r2];break;}
				case 't': {mp[j][r1]=triplet[3][r2];break;}
				default:
					{
						printf("wrong letter %c in %d position\n",mp[j][r1],j);mp[j][r1]='n';
					}
				}		
		}
	}	
}
void town::get_copy(town *a, int *n_mut, int n_olig)
{
	int i,j;
	a->bpv_n=bpv_n;
	a->bpv_wy=bpv_wy;
	a->bpv_wn=bpv_wn;
	a->bmat=bmat;			
	a->bpos=bpos;			
	a->ori=ori;		
	a->fit=fit;
	strcpy(a->seq,seq);
	for(j=0;j<n_olig;j++)
	{
		for(i=0;i<n_mut[j];i++)a->mp[j][i]=mp[j][i];
	}		
}
void town::print_all(int *n_mut, int n_olig, FILE* out)
{
 int i,j;
 fprintf(out,"F %.6f PwY %.6f PwN %.6f PN %.6f M %3d P %3d%c ",fit,bpv_wy,bpv_wn,bpv_n,bmat,bpos,ori);
for(j=0;j<n_olig;j++)
{
	fprintf(out, " ");
	for(i=0;i<n_mut[j];i++)fprintf(out, "%c",mp[j][i]);
}
fprintf(out, " %s\n",seq);
}
int compare_pop( const void *X1, const void *X2 )
{
	struct town *S1 = (struct town *)X1;
	struct town *S2 = (struct town *)X2; 
	if(S1->fit - S2->fit >0)return 1;
	if(S1->fit - S2->fit <0)return -1;		
	return 0;
}
int compare_qq( const void *X1, const void *X2 )//otbor men'6ego
{
	double X=(*(double*)X1 - *(double*)X2 );
	if(X>0)return -1;
	if(X<0)return 1;		
	return 0;
}
struct matr {
	int len;
	double min;//min sco
	double raz;// (max-min) sco
	double **wes;
} *pwm;
struct pval {
	int nthr;
	double *thr;
	double *lgpv;
} *pvt;
int GomTown(town a, town b, int *n_mut, int n_olig)
{	
	int i,j;	
	for(i=0;i<n_olig;i++)
	{
		for(j=0;j<n_mut[i];j++)if(b.mp[i][j]!=a.mp[i][j])return 0;
	}	
	return 1;
}
int ComplStr(char *d)
{
	char d1[LENGTH];
	int i, len;
	len=strlen(d);
	strcpy(d1,d);
//	memset(d,0,sizeof(d));
	for(i=0;i<len;i++)
	{
		switch(d1[len-i-1])
		{
			case 'a':{d[i]='t';break;}
			case 't':{d[i]='a';break;}
			case 'c':{d[i]='g';break;}
			case 'g':{d[i]='c';break;}
			case 'A':{d[i]='T';break;}
			case 'T':{d[i]='A';break;}
			case 'C':{d[i]='G';break;}
			case 'G':{d[i]='C';break;}
			case 'N':{d[i]='N';break;}
			case 'n':{d[i]='n';break;}
			default: d[i]='n';
		}
	}
	return 1;
}
double EvalFit(town *a, int **mutpos, char *olig, char *olig0, int *n_mut, int n_olig, char *mname, int n_pwm, char *mpv_ext, int m_positive, int n_pos_sites, int *m_ignore)
{	
	int i,j,m,k,p, ori;
	double wfit[25]={1,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125,0.00390625,0.001953125,0.000976563,0.000488281,0.000244141,	0.00012207,
		6.10352E-05,3.05176E-05,1.52588E-05,7.62939E-06,3.8147E-06,1.90735E-06,9.53674E-07,4.76837E-07,2.38419E-07,1.19209E-07,5.96046E-08};	
	double sum_rfit=1.99999994;
	double rfit[NMATRICES];	
	double mpos_score[2*LENGTH];
	int n_pos_cases=0;	
	char str[2][LENGTH], str0[2][LENGTH];
	
	for(m=0;m<n_pwm;m++)rfit[m]=0;
	strcpy(str[0],olig);
	strcpy(str0[0],olig0);
	int len=strlen(str[0]);
	for(i=0;i<n_olig;i++)
	{
		for(j=0;j<n_mut[i];j++)str[0][mutpos[i][j]]=a->mp[i][j];
	}	
	strcpy(str[1],str[0]);
	strcpy(str0[1],str0[0]);
	ComplStr(str[1]);
	ComplStr(str0[1]);
	double ret=0;	
	for(m=0;m<n_pwm;m++)
	{
		int ignore=0;
		for(k=0;m_ignore[k]!=-1;k++)
		{
			if(m==m_ignore[k])
			{
				ignore=1;
				break;
			}
		}
		if(ignore==1)continue;
		//printf("%d ",m+1);
		/*if(m==m_positive)
		{
			rfit[m]=0;
			continue;
		}*/
		int lenp=len-pwm[m].len+1;
		double best_sco=0;
		int best_pos;
		char best_ori, best_seq[PWMLEN];
		int last_num = pvt[m].nthr-1;
		double thresh=pvt[m].thr[last_num];
		for(ori=0;ori<2;ori++)
		{
			for(k=0;k<lenp;k++)
			{
				double sco=0;
				for(p=0;p<pwm[m].len;p++)
				{
					switch (str[ori][k+p]){
					case 'a':{sco+=pwm[m].wes[p][0];break;}
					case 'c':{sco+=pwm[m].wes[p][1];break;}
					case 'g':{sco+=pwm[m].wes[p][2];break;}
					case 't':{sco+=pwm[m].wes[p][3];break;}
					default: {break;}
					}
				}
				sco-=pwm[m].min;
				sco/=pwm[m].raz;
				if(m==m_positive)
				{
					mpos_score[n_pos_cases++]=sco;					
				}
				if(sco<=thresh)continue;
				if(sco>best_sco)
				{
					best_sco=sco;					
					strncpy(best_seq,&str0[ori][k],pwm[m].len);
					best_seq[pwm[m].len]='\0';
					if(ori==0)
					{
						best_ori='+';
						best_pos=k+1;
					}
					else 
					{
						best_ori='-';
						best_pos=lenp-k-1;
					}
				}				
			}
		}		
		if(best_sco < thresh)
		{
			/*if (m == m_positive)
			{
				a->fit = 10;
				return a->fit;
			}*/
			rfit[m]=0;		
			continue;
		}
		if(m!=m_positive)
		{
			if(best_sco >= pvt[m].thr[0])
			{						
				rfit[m]=pvt[m].lgpv[0];
				if(ret<rfit[m])
				{
					ret=rfit[m];
					a->bpos=best_pos;
					a->bmat=m;
					a->ori=best_ori;
					strncpy(a->seq,best_seq,pwm[m].len);
					a->seq[pwm[m].len]='\0';
					continue;
				}
			}
			else
			{								
				for(k=last_num;k>=1;k--)
				{
					if(best_sco>=pvt[m].thr[k] && best_sco<pvt[m].thr[k-1])
					{
						rfit[m]=pvt[m].lgpv[k];																				
						if(ret<rfit[m])
						{
							ret=rfit[m];
							a->bpos=best_pos;
							a->bmat=m;
							a->ori=best_ori;
							strncpy(a->seq,best_seq,pwm[m].len);
							a->seq[pwm[m].len]='\0';
							break;
						}
					}
				}				
			}
		}
		else rfit[m]=0;
	}		
	if((m_positive>=0 && m_positive<n_pwm) && n_pos_cases<n_pos_sites)return 0;
	a->bpv_n=ret;	
	/*
	for(m=0;m<n_pwm;m++)
	{
		printf("%.3f\t",rfit[m]);
		if((m+1)%10==0)printf("\n");
	}*/
	qsort(rfit,n_pwm,sizeof(double),compare_qq);	
	/*for(m=0;m<n_pwm;m++)
	{
		printf("%.3f\t",rfit[m]);
		if((m+1)%10==0)printf("\n");
	}*/
	double wei_n=0, wei_y=0;
	for(i=0;i<25;i++)wei_n+=rfit[i]*wfit[i];
	wei_n/=sum_rfit;
	a->bpv_wn=wei_n;
	if(m_positive<0 || m_positive>=n_pwm)a->bpv_wy=1;
	else
	{
		qsort(mpos_score,n_pos_cases,sizeof(double),compare_qq);	
		for(i=0;i<n_pos_sites;i++)wei_y+=mpos_score[i];
		wei_y/=n_pos_sites;
		if(wei_y >= pvt[m_positive].thr[0])a->bpv_wy=pvt[m_positive].lgpv[0];
		else
		{
			int last_num = pvt[m_positive].nthr-1;
			if(wei_y < pvt[m_positive].thr[last_num])a->bpv_wy=pvt[m_positive].lgpv[last_num];
			else
			{
				for(k=last_num;k>=1;k--)
				{
					if(wei_y>=pvt[m_positive].thr[k] && wei_y<pvt[m_positive].thr[k-1])
					{
						a->bpv_wy=pvt[m_positive].lgpv[k];
						break;
					}
				}
			}
		}
	}
	a->fit=a->bpv_wn/a->bpv_wy;	
	return a->fit;
}
int MutPoint(town *a, int *n_mut, int n_olig, int **mut_pos)
{
	int i, j;
	if(a->bmat==0){j=rand()%n_olig;return j;}
	else
	{
		int s1=a->bpos, s2=s1+pwm[a->bmat].len-1;
		int jmut[NOLIG], nj=0;
		for(i=0;i<n_olig;i++)jmut[i]=-1;
		int ii=0;		
		for(i=0;i<n_olig;i++)
		{
			int n_mut1=n_mut[i]-1;
			int m1 =mut_pos[i][0], m2=mut_pos[i][n_mut1];
			if(m1<=s1 && s1 <=m2){jmut[ii++]=i;continue;}
			if(m1<=s2 && s2 <=m2){jmut[ii++]=i;continue;}
			if(s1<=m1 && m1 <=s2){jmut[ii++]=i;continue;}
			if(m1<=s2 && s2 <=m2){jmut[ii++]=i;continue;}
		}
		if(ii>0)
		{
			int r3=rand()%ii;
			j=jmut[r3];
			return j;
		}
	}	
	return -1;
}
int Mut1(town *a, int *n_mut, int n_olig, char **wt_nt, int j, int *n_mut_deg)
{
	int i,k,s;
	char triplet[4][4]={"cgt","agt","act","acg"};	// a c g t
	for(i=0;i<n_mut[j];i++)a->mp[j][i]=wt_nt[j][i];	
	char muts[LENGTH];
	int n_mut_deg_cur=0;
	for(k=0;k<n_mut_deg[j];k++)
	{
		int gom=1, r1;
		do
		{
			gom=1;
			r1=rand()%n_mut[j];			
			for(s=0;s<n_mut_deg_cur;s++)
			{
				if(muts[s]==r1)
				{
					gom=0;
					break;
				}
			}
		}
		while(gom==0);
		muts[n_mut_deg_cur]=r1;
		n_mut_deg_cur++;
		int r2=rand()%3;
		switch (a->mp[j][r1]) {
			case 'a': {a->mp[j][r1]=triplet[0][r2];break;}
			case 'c': {a->mp[j][r1]=triplet[1][r2];break;}
			case 'g': {a->mp[j][r1]=triplet[2][r2];break;}
			case 't': {a->mp[j][r1]=triplet[3][r2];break;}
			default:{printf("wrong letter %c in %d position\n",a->mp[j][r1],j);return -1;}
			}		
	}
	return 1;
}
int Reco2(town *a, town *b, int *n_mut, int n_olig, int j)
{
	int i;	
	int gom=1;		
	for(i=0;i<n_mut[j];i++)
	{		
		if(a->mp[j][i]!=b->mp[j][i])
		{
			gom=0;
			break;
		}
	}
	if(gom==1)
	{
		return -1;	
	}
	else
	{
		for(i=0;i<n_mut[j];i++)
		{
			char x = a->mp[j][i];			
			a->mp[j][i] = b->mp[j][i];
			b->mp[j][i] = x;			
		}
		return 1;		
	}	
}
int StrNStr(char *str,char c, int n)
{
	if(n==0)return -1;
	int i, len=strlen(str);
	int k=1;
	for(i=0;i<len;i++)
	{
		if(str[i]==c)
		{
			if(k==n)return i;
			k++;
		}
	}
	return -1;
}
int UnderStol(char* str, int nstol, char* ret, size_t size, char sep)
{
	memset(ret, 0, size);
	int p1, p2, len;
	if (nstol == 0)
	{
		p2 = StrNStr(str, sep, 1);
		if (p2 == -1)p2 = strlen(str);
		strncpy(ret, str, p2);
		ret[p2] = '\0';
		return 1;
	}
	else
	{
		p1 = StrNStr(str, sep, nstol);
		p2 = StrNStr(str, sep, nstol + 1);
		if (p2 == -1)
		{
			p2 = strlen(str);
		}
		if (p1 == -1 || p2 == -1) return -1;
		len = p2 - p1 - 1;
		strncpy(ret, &str[p1 + 1], len);
		ret[len] = '\0';
		return 1;
	}
}
void DelChar(char *str,char c)
{
	int i, lens, size;

	size=0;
	lens=strlen(str);
	for(i=0;i<lens;i++)
	{
		if(str[i]!=c)str[size++]=str[i];
	}
	str[size]='\0';
}
char *TransStr(char *d)
{
   int i, c, lens;
   lens=strlen(d);
   for(i=0;i<lens;i++)
   {
	c=int(d[i]);
	if(c<97) d[i]=char(c+32);
	//else break;
   }
   return(d);
}
char TransStrBack(char &d)//a->A
{
   int c=int(d);
	if(c>=97) d=char(c-32);
   return(d);
}
int main(int argc, char *argv[])
{
	int i,j,k,m;
	char olig[LENGTH], olig0[LENGTH], str[1000], filei_seq[300], filei_tab[300], filei_mat[300], fileo[300], file_log[300], path_pos[300], path_neg[300];
	char mname[50], mext[]=".pwm", mpv_ext[]=".dist";
	FILE *in, *out, *outlog;
	if(argc!=11)
	{
		puts("Syntax: 1path pos_matrix 2path neg_matrices 3file input_sequence 4file_razmetka_mutaciy 5int matrix_count\n");
		puts("6char matrix_name 7file out 8int anchor_mode (1,0 yes,no) 9double probability_mutation(0..1) 10file log");
		return -1;
	}
	strcpy(path_pos, argv[1]);
	strcpy(path_neg, argv[2]);
	strcpy(filei_seq, argv[3]);
	strcpy(filei_tab, argv[4]);
	int n_pwm=atoi(argv[5]);
	strcpy(mname,argv[6]);
	strcpy(fileo,argv[7]);
	int anchor=atoi(argv[8]);
	double prob_mut=atof(argv[9]);
	strcpy(file_log, argv[10]);
	int m_positive;
	int m_ignore_count = 0;
	//int m_ignore[] = { 2, 82, 83, 84, 102, 183, 184, 186, 207, 212, 213, 217, 286, 299, 300, 303, 439, -1 };// 2020 from mcot, 82,83,84 = ARFs
	//int m_ignore[] = { 2, 34, 35, 36, 37, 41, 75, 81, 102, 183, 184, 186, 207, 212, 213, 217, 286, 299, 300, 303, 439, -1 };//2021 cbf3
//int m_ignore[]={ 2, 102, 183, 184, 186, 207, 212, 213, 217, 227, 286, 299, 300, 303, 439, -1 };// 2020 from mcot, 227 = EIN3		
	//int m_ignore[] = { 2, 4, 102, 183, 184, 186, 207, 212, 213, 217, 227, 286, 299, 300, 303, 439, -1 };// 2023, 4 = FUS3, 227 = EIN3
	//int m_ignore[] = { 2, 102, 183, 184, 186, 207, 212, 213, 214, 217, 229, 286, 299, 300, 303, 439, -1 };//2021 x = 214 camta1,229 far
	//int n_pos_sites=20;//arf number of positive sites to take into account
	int n_pos_sites = 10;//ein3
	int m_ignore[] = { 2, 102, 183, 184, 186, 207, 212, 213, 217, 227, 286, 299, 300, 303, 439, -1 };// no anchor 2024
	if(anchor==1)
	{
		m_positive=0;
		for (i = 0; m_ignore[i] != -1; i++)m_ignore_count++;
	//	m_ignore[0]=227, m_ignore[1]=228,  m_ignore[2]=-1;//ein3
	//	m_ignore[0]=82, m_ignore[1]=83,m_ignore[2]=84,m_ignore[3]=-1; //arf
	}
	else
	{
		m_positive=-1;
		m_ignore[0] = 1000;// , m_ignore[1] = -1;
	}
	if((in=fopen(filei_seq,"rt"))==NULL)
	{
		 printf("Input file %s can't be opened!",filei_seq);
		 return -1;
	}
	fgets(olig,sizeof(olig),in);//head
	fgets(olig,sizeof(olig),in);
	DelChar(olig,'\n');
	DelChar(str, '\r');
	strcpy(olig0,olig);
	TransStr(olig);
	fclose(in);
	if((in=fopen(filei_tab,"rt"))==NULL)
	{
		 printf("Input file %s can't be opened!",filei_tab);
		 return -1;
	}
	int n_olig=0;//  4islo blokov mutirovaniya
	while(fgets(str,sizeof(str),in)!=NULL)
	{
		if(isdigit(str[0]))n_olig++;
	}
	rewind(in);
	int *n_mut;	
	n_mut = new int[n_olig];
	if(n_mut==NULL){puts("Out of memory...");return-1;}				
	int *n_mut_deg;	
	n_mut_deg = new int[n_olig];
	if(n_mut_deg==NULL){puts("Out of memory...");return-1;}				
	for(j=0;j<n_olig;j++)
	{
		if(fgets(str,sizeof(str),in)!=NULL)
		{
			n_mut[j]=0;
			if(isdigit(str[0]))
			{
				n_mut[j]++;
			}
			int len=strlen(str);
			char tab='\t';
			for(i=0;i<len;i++)
			{
				if(str[i]==tab)n_mut[j]++;
			}
		}
	}
	rewind(in);
	{
		int n_mut_max=0;
		for(j=0;j<n_olig;j++)
		{
			if(n_mut[j]>n_mut_max)n_mut_max=n_mut[j];
		}
		if(n_olig>NOLIG || n_mut_max>NMUT)
		{
			printf("Wrong parameters %d no.of oligs, %d no.of mutations\n",n_olig,n_mut_max);
			return-1;
		}
	}
	for(j=0;j<n_olig;j++)
	{
		int n_mut_deg1=(int)(prob_mut*n_mut[j]);
		n_mut_deg[j]=Max(1,n_mut_deg1);
	}
	int n_mut_total=0;
	for(j=0;j<n_olig;j++)n_mut_total+=n_mut[j];
	//int n_mut_enter=(int)(n_mut_total*prob_mut);
	int **mutpos;
	mutpos = new int*[n_olig];
	if(mutpos==NULL){puts("Out of memory...");return-1;}			
	for(k=0;k<n_olig;k++)
	{
		mutpos[k] = new int[n_mut[k]];
		if(mutpos[k]==NULL ){puts("Out of memory...");return-1;}
	}
	char **wt_nt;
	wt_nt = new char*[n_olig];
	if(wt_nt==NULL){puts("Out of memory...");return-1;}			
	for(k=0;k<n_olig;k++)
	{
		wt_nt[k] = new char[n_mut[k]];
		if(wt_nt[k]==NULL ){puts("Out of memory...");return -1;}
	}
	for(j=0;j<n_olig;j++)
	{
		if(fgets(str,sizeof(str),in)==NULL)
		{
			printf("Unexpected end of %s file\n",filei_tab);
			return -1;
		}
		char value[10], razd='\t';
		int sta=0, i_num=0, iend=strlen(str);	
		for(i=0;i<=iend;i++)
		{
			if(i==iend || str[i]==razd)
			{
				value[sta]='\0';
				mutpos[j][i_num]=atoi(value);
				mutpos[j][i_num++]--;
				memset(value,0,sizeof(value));
				sta=0;
			}
			else value[sta++]=str[i];
		}
	}
	fclose(in);
	for(i=0;i<n_olig;i++)
	{
		for(j=0;j<n_mut[i];j++)wt_nt[i][j]=olig[mutpos[i][j]];
	}
	char sep = '\t';
//matricy
	pwm = new matr[n_pwm];
	if(pwm==NULL ){puts("Out of memory...");return -1;}
	for(m=0;m<n_pwm;m++)
	{
		if (m == 0)strcpy(filei_mat, path_pos);
		else strcpy(filei_mat, path_neg);
		strcat(filei_mat, mname);		
		char bufn[10];
		memset(bufn,0,sizeof(bufn));
		sprintf(bufn, "%d", m);
		strcat(filei_mat,bufn);
		strcat(filei_mat,mext);
		if((in=fopen(filei_mat,"rt"))==NULL)
		{
			printf("Input file %s can't be opened!",filei_mat);
			return -1;
		} 
		pwm[m].len=0;
		while(fgets(str,sizeof(str),in)!=NULL)
		{	
			if(isdigit(str[0]) || str[0]=='-')pwm[m].len++;
		}
		rewind(in);
		pwm[m].wes = new double*[pwm[m].len];
		if(pwm[m].wes==NULL){puts("Out of memory...");return -1;}			
		for(k=0;k<pwm[m].len;k++)
		{
			pwm[m].wes[k] = new double[4];
			if(pwm[m].wes[k]==NULL ){puts("Out of memory...");return -1;}
		}
		pwm[m].raz=pwm[m].min=0;
		fgets(str,sizeof(str),in);//head
		for(i=0;i<pwm[m].len;i++)		
		{
			if(fgets(str,sizeof(str),in)==NULL)			
			{
				printf("Unexpected end of %s file\n",filei_mat);
				return -1;
			}
			if(isdigit(str[0]) || str[0]=='-')
			{
				char sco[20];
				double pwmmin=100;
				double pwmmax=-100;
				for(j=0;j<4;j++)
				{
					int u = UnderStol(str, j, sco, sizeof(sco), sep);
					if(u==-1)
					{				
						printf("Error reading %s file\n",filei_mat);
						return -1;
					}
					double val = atof(sco);		
					if(val<pwmmin)pwmmin=val;
					if(val>pwmmax)pwmmax=val;
					pwm[m].wes[i][j]=val;
					//printf("%f\t",wei[i][j]);
				}
				pwm[m].raz+=pwmmax;
				pwm[m].min+=pwmmin;
				//printf("\n");
			}
		}
		pwm[m].raz-=pwm[m].min;
		fclose(in);
	}
	//porogi
	pvt = new pval[n_pwm];
	if(pvt==NULL ){puts("Out of memory...");return -1;}
	for(m=0;m<n_pwm;m++)
	{
		if (m == 0)strcpy(filei_mat, path_pos);
		else strcpy(filei_mat, path_neg);
		strcat(filei_mat, mname);		
		char bufn[10];
		memset(bufn,0,sizeof(bufn));
		sprintf(bufn, "%d", m);
		strcat(filei_mat,bufn);
		strcat(filei_mat,mpv_ext);
		if((in=fopen(filei_mat,"rt"))==NULL)
		{
			printf("Input file %s can't be opened!",filei_mat);
			return -1;
		} 		
		int n_str=0;
		double pv_thr=3;
		while(fgets(str,sizeof(str),in)!=NULL)
		{	
			char sco[20];
			if(isdigit(str[0]))
			{				
				int u = UnderStol(str, 1, sco, sizeof(sco), sep);
				if(u==-1)
				{				
					printf("Error reading %s file\n",filei_mat);
					return -1;
				}
				double pv=atof(sco);
				if(pv>pv_thr)n_str++;
				else {n_str++;break;}
			}
		}
		pvt[m].nthr=Min(n_str,NPVT);
		pvt[m].lgpv = new double[pvt[m].nthr];
		if(pvt[m].lgpv==NULL){puts("Out of memory...");return -1;}			
		pvt[m].thr = new double[pvt[m].nthr];
		if(pvt[m].thr==NULL){puts("Out of memory...");return -1;}		
		int str_gap;
		if(n_str<=NPVT)str_gap=0;
		else str_gap=n_str/NPVT;
		rewind(in);
		for(i=0;i<pvt[m].nthr;i++)		
		{
			if(fgets(str,sizeof(str),in)==NULL)			
			{
				printf("Unexpected end of %s file\n",filei_mat);
				return -1;
			}
			if(isdigit(str[0]))
			{
				char sco[20];				
				int u = UnderStol(str, 1, sco, sizeof(sco), sep);
				if(u==-1)
				{				
					printf("Error reading %s file\n",filei_mat);
					return -1;
				}
				double val=atof(sco);			
				pvt[m].lgpv[i] = val;// -log10(val);
				pvt[m].thr[i]=atof(str);				
			}
			for(k=0;k<str_gap;k++)
			{
				if(fgets(str,sizeof(str),in)==NULL)			
				{
					printf("Unexpected end of %s file\n",filei_mat);
					return -1;
				}
			}
		}
		fclose(in);
		/*if((m+1)%200==0)
		{
			int yy=0;
		}*/
		//printf("Matrix %d Nthr %5d Thr %.6f Pvalue %g\n",m+1,pvt[m].nthr,pvt[m].thr[pvt[m].nthr-1],pow((double)10,-pvt[m].lgpv[pvt[m].nthr-1]));
	}	
	if ((outlog = fopen(file_log, "wt")) == NULL)
	{
		fprintf(outlog, "Input file %s can't be opened!\n", file_log);
		exit(1);
	}
	dmut.init0(wt_nt,n_mut,n_olig);
	EvalFit(&dmut,mutpos,olig,olig0,n_mut,n_olig,mname,n_pwm,mpv_ext,m_positive,n_pos_sites,m_ignore);
	dmut.print_all(n_mut,n_olig,outlog);
	srand( (unsigned)time( NULL ) );
	int rec_win =2;
	int m_success_no_max=5, r_success_no_max=1, r_cross_no_max=10, total_m_success, total_r_success;
	for(i=0;i<POOL;i++)
	{		
	//	printf("\b\b\b%3d",i);
		int gom=0;
		do
		{
			dmut.init(wt_nt,n_mut,n_olig,n_mut_deg);
			//dmut.print_all(n_mut,n_olig);
			for(j=0;j<i;j++)
			{
				if(GomTown(dmut,pop[j],n_mut,n_olig)==1)
				{
					gom=1;
					break;
				}
			}
		}
		while(gom==1);
		//dmut.print_all(n_mut,n_olig);
		dmut.get_copy(&pop[i],n_mut,n_olig);		
		EvalFit(&pop[i],mutpos,olig,olig0,n_mut,n_olig,mname,n_pwm,mpv_ext,m_positive,n_pos_sites,m_ignore);		
	}
	qsort((void*)pop,POOL, sizeof(pop[0]), compare_pop);	
	for(i=0;i<ELIT;i++)pop[i].print_all(n_mut,n_olig,outlog);
	int iter=0;
	double fit0=pop[0].fit;
	int empty_cycle=0;
	do
	{
		total_m_success=0;
		total_r_success=0;
		fprintf(outlog, "Mutation %d    \n",iter+1);
		for(i=0;i<POOL;i++)
		{
			int m_success_no=0;
			//if((i+1)%10==0)printf("\b\b\b%3d",i+1);
			if(pop[i].fit==0)continue;
			do
			{
				pop[i].get_copy(&dmut,n_mut,n_olig);
//				dmut.print_all(n_mut, n_olig);
				int x=MutPoint(&dmut,n_mut,n_olig,mutpos);
			//	int x=rand()%n_olig;
				int mut = -1;
				if(x!=-1)Mut1(&dmut, n_mut, n_olig, wt_nt, x, n_mut_deg);
	//			dmut.print_all(n_mut, n_olig);
				if(mut==-1)
				{
					m_success_no++;
					continue;
				}
				int gom=0;
				for(j=0;j<POOL;j++)
				{
					if(GomTown(dmut,pop[j],n_mut,n_olig)==1){gom=1;break;}
				}
				if(gom==1)
				{
					m_success_no++;
					continue;
				}
				//	dmut.print_all(size,max);				
				EvalFit(&dmut,mutpos,olig,olig0,n_mut,n_olig,mname,n_pwm,mpv_ext,m_positive,n_pos_sites,m_ignore);
				if(dmut.fit<pop[i].fit)
				{
					dmut.get_copy(&pop[i],n_mut,n_olig);	
					m_success_no=0;
					total_m_success++;
				}
				else m_success_no++;
			}
			while(m_success_no<m_success_no_max);
		}
		qsort((void*)pop,POOL, sizeof(pop[0]), compare_pop);	
		fprintf(outlog, "GA %d\tMut %d\n",iter+1,total_m_success);
		for(i=0;i<ELIT;i++)pop[i].print_all(n_mut,n_olig,outlog);	
		fprintf(outlog,"Recombination %d      ",iter+1);
		for(i=0;i<POOL;i++)
		{			
			//if((i+1)%10==0)			
		//	if ((i + 1) % 20 == 0)printf("\b\b\b\b\b\b\b%4d", i);
			if(pop[i].fit==0)continue;
	//		printf("Rec 1st %3d\t",i);			
			for(j=i+1;j<POOL;j++)
			{				
			//	if(i==j)continue;
				if(pop[j].fit==0)continue;
		//		printf("Rec 2nd %3d\t",j);						
				int r_success_no=0;
				do
				{
					int rec;
					pop[i].get_copy(&drec[0],n_mut,n_olig);			
					pop[j].get_copy(&drec[1],n_mut,n_olig);	
				//	printf("\n");
				//	for(k=0;k<2;k++)drec[k].print_all(n_mut,n_olig);
					int x, det;
					if(drec[0].fit>drec[1].fit)det=0;
					else
					{
						if(drec[0].fit<drec[1].fit)det=1;
						else det=rand()%2;
					}
					x=MutPoint(&drec[det],n_mut,n_olig,mutpos);
					rec = Reco2(&drec[0],&drec[1],n_mut,n_olig,x);
					//for(k=0;k<2;k++)drec[k].print_all(n_mut,n_olig);
					if(rec==-1)
					{
						r_success_no++;
						continue;
					}
					int gom=0;
					for(k=0;k<2;k++)
					{						
						for(m=0;m<POOL;m++)
						{
							if(GomTown(drec[k],pop[m],n_mut,n_olig)==1){gom=1;break;}
						}
						if(gom==1)break;
					}
					if(gom==1)
					{
						r_success_no++;
						break;
					}
					double fitmin=Min(pop[i].fit,pop[j].fit);
					for(k=0;k<2;k++)EvalFit(&drec[k],mutpos,olig,olig0,n_mut,n_olig,mname,n_pwm,mpv_ext,m_positive,n_pos_sites,m_ignore);
					if(drec[0].fit<fitmin || drec[1].fit<fitmin)
					{
					//	for(k=0;k<2;k++)drec[k].print_all(size,max);	
						drec[0].get_copy(&pop[i],n_mut,n_olig);			
						drec[1].get_copy(&pop[j],n_mut,n_olig);	
						total_r_success++;
					}
					else r_success_no++;										
				}
				while(r_success_no<r_success_no_max);
			}
		}		
		fprintf(outlog,"\n");
		qsort((void*)pop,POOL, sizeof(pop[0]), compare_pop);
		for(i=0;i<ELIT;i++)pop[i].print_all(n_mut,n_olig,outlog);
		iter++;
		double fit0_here=pop[0].fit;
		if(fit0==fit0_here)empty_cycle++;
		else empty_cycle=0;
		fit0=fit0_here;
		fprintf(outlog, "GA %d\tMut %d\tRec %d\tEmpty_cycle %d\n",iter+1,total_m_success,total_r_success,empty_cycle);
		fclose(outlog);
		if ((outlog = fopen(file_log, "at")) == NULL)
		{
			fprintf(outlog, "Input file %s can't be opened!\n", file_log);
			exit(1);
		}
	}
	while(total_m_success+total_r_success>POOL && empty_cycle<10);
	fprintf(outlog, "Final\n");
	fclose(outlog);
	//int ext_elit= Min(POOL, 5*ELIT);
	for(i=0;i<ELIT;i++)pop[i].print_all(n_mut,n_olig,outlog);
	if((out=fopen(fileo,"wt"))==NULL)
	{
		printf("Output file %s can't be opened!\n", fileo);
		return -1;
	}	
	fprintf(out,">WT variant\n%s\n",olig0);
	for(k=0;k<POOL;k++)
	{
		strcpy(str,olig);
		fprintf(out,">%d\t",k+1);
		fprintf(out,"Fit %.6f\tPvNeg %.6f\tPvNegW %.6f\tPvPosW %.6f\t",pop[k].fit,pop[k].bpv_n,pop[k].bpv_wn, pop[k].bpv_wy);
		fprintf(out,"M %d\tP %d%c %s",pop[k].bmat,pop[k].bpos,pop[k].ori,pop[k].seq);
		for(j=0;j<n_olig;j++)
		{			
			for(i=0;i<n_mut[j];i++)
			{			
				char c=pop[k].mp[j][i];
				if(c!=wt_nt[j][i])TransStrBack(c);
				str[mutpos[j][i]]=c;
			}
		}
		fprintf(out,"\n%s\n",str);
	}
   strcpy(str,olig);
   fclose(out);
}