#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define NOLIG 1050 // max number of oligonucleotides 
#define PWMLEN 32//max length of PWM
#define LENGTH 350 //dlina polimera = 21*10 
#define NPVT 7000 // 4islo porogov zadannyh p-value
#define POOL 400 // size of population = 4islo osobey populyacii
#define ELIT 25 // size of population = 4islo osobey populyacii - na vyhod elita
#define NMATRICES 600

char **seq;
char **seq0;

void Mix(int *a, int *b)// bez commentov
{
  int buf=*a;
  *a=*b;
  *b=buf;
}
void BigMix1(int *d1, int len) // pereme6ivanie stroki
{
    int r;    
 for(r=0;r<len-1;r++)
 {
	 Mix(&d1[r], &d1[1+r+(rand()%(len-1-r))]);	 
 }
}

struct town {
	int da[NOLIG];// eti est' = order
	int deg[NOLIG];// zaselennost
	double fit; //fitness
	double bpv_n;// best score negative
	double bpv_wn;//best score negative weighted
	double bpv_wy;// av best score positive
	int bpos;
	int bmat;
	char ori;
	char seq[PWMLEN];
	void init(int size, int max);
	void get_copy(town *a ,int size, int max);
	void print_all(int size, int max, FILE *out);	
	int check(int size, int max);
} dmut, drec[2], pop[POOL];
void town::init(int size, int max)
{
	int i;
	fit=0;
	bpv_wn = bpv_wy = bpv_n = 0;
    bmat=0, bpos=0, ori='+';
	memset(seq,'\0',sizeof(seq));
	int arr[NOLIG], inx[NOLIG];
	for(i=0;i<size;i++)arr[i]=1;
	for(i=size;i<max;i++)arr[i]=0;
	BigMix1(arr,max);		
	int i1=0; //zanyat nomer
	for(i=0;i<max;i++)
	{
		if(arr[i]==1)inx[i1++]=i;
	}
	BigMix1(inx,size);
	for(i=0;i<size;i++)da[i]=inx[i];	
	fit=0;
	for(i=0;i<max;i++)deg[i]=0;
	for(i=0;i<size;i++)deg[da[i]]=1;	
}
int town::check(int size, int max)
{
	int i,j;
	for(i=0;i<size;i++)
	{		
		if(deg[da[i]]!=1)
		{
			return -1;
		}
	}
	for(i=0;i<max;i++)
	{
		if(deg[i]==0)
		{
			for(j=0;j<size;j++)
			{
				if(da[j]==i)return -1;
			}
		}
	}	
	for(i=0;i<size;i++)
	{
		int dai=da[i];
		for(j=i+1;j<size;j++)
		{
			if(dai==da[j])return -1;
		}
	}
	int sum=0;
	for(i=0;i<max;i++)sum+=deg[i];
	if(sum!=size)return -1;
	return 1;
}
void town::get_copy(town *a ,int size, int max)
{
	int i;
	a->bmat=bmat;			
	a->bpv_n = bpv_n;
	a->bpv_wy = bpv_wy;
	a->bpv_wn = bpv_wn;
	a->bpos=bpos;			
	a->ori=ori;		
	a->fit=fit;			
	strcpy(a->seq,seq);
	for(i=0;i<size;i++)a->da[i]=da[i];
	for(i=0;i<max;i++)a->deg[i]=deg[i];			
}
void town::print_all(int size, int max, FILE* out)
{
 int i; 
 fprintf(out,"F %.6f PwY %.6f PwN %.6f PN %.6f M %3d P %3d%c ", fit, bpv_wy, bpv_wn, bpv_n, bmat, bpos, ori);
 for(i=0;i<size;i++)fprintf(out,"%2d ",da[i]+1);
/* printf(" deg: ");
 for(i=0;i<max;i++)
 {	
	printf("%1d ",deg[i]);
	if((i+1)%5==0)printf(" ");
 }*/
fprintf(out," %s\n",seq);
}
int compare_pop( const void *X1, const void *X2 )
{
	struct town *S1 = (struct town *)X1;
	struct town *S2 = (struct town *)X2; 
	if(S1->fit - S2->fit >0)return 1;
	if(S1->fit - S2->fit <0)return -1;		
	return 0;
}
int compare_qq( const void *X1, const void *X2 )
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
int GomTown(town a, town b, int size)
{	
	int i;	
	for(i=0;i<size;i++)
	{
		if(b.da[i]!=a.da[i])return 0;
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
//double EvalFit(town *a, int **mutpos, char *olig, char *olig0, int *n_mut, int n_olig, char *mname, int n_pwm, char *mpv_ext, int m_positive, int n_pos_sites, int *m_ignore)
double EvalFit(town *a, int olen, int size, char *mname, int n_pwm, char *mpv_ext, int m_positive, int n_pos_sites, int *m_ignore)
{
	int i, m, k, p, ori;
	char str[2][LENGTH], str0[2][LENGTH];
	double wfit[25] = { 1,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125,0.00390625,0.001953125,0.000976563,0.000488281,0.000244141,	0.00012207,
		6.10352E-05,3.05176E-05,1.52588E-05,7.62939E-06,3.8147E-06,1.90735E-06,9.53674E-07,4.76837E-07,2.38419E-07,1.19209E-07,5.96046E-08 };
	double mpos_score[2 * LENGTH];
	int n_pos_cases = 0;
	double sum_rfit = 1.99999994;

	double rfit[NMATRICES];
	for (m = 0; m < n_pwm; m++)rfit[m] = 0;

	memset(str[0], '\0', sizeof(str[0]));
	memset(str[1], '\0', sizeof(str[1]));
	for (i = 0; i < size; i++)strcat(str[0], seq[a->da[i]]);
	int len = olen * size;
	strcpy(str[1], str[0]);
	ComplStr(str[1]);
	str[1][len] = str[0][len] = '\0';

	memset(str0[0], '\0', sizeof(str0[0]));
	memset(str0[1], '\0', sizeof(str0[1]));
	for (i = 0; i < size; i++)strcat(str0[0], seq0[a->da[i]]);
	strcpy(str0[1], str0[0]);
	ComplStr(str0[1]);
	str0[1][len] = str0[0][len] = '\0';

	double ret = 0;	

	for (m = 0; m < n_pwm; m++)
	{
		int ignore = 0;
		for (k = 0; m_ignore[k] != -1; k++)
		{
			if (m == m_ignore[k])
			{
				ignore = 1;
				break;
			}
		}
		if (ignore == 1)continue;
		//printf("%d ",m);
		/*if(m==m_positive)
		{
			rfit[m]=0;
			continue;
		}*/
		int lenp = len - pwm[m].len + 1;
		double best_sco = 0;
		int best_pos;
		char best_ori, best_seq[PWMLEN];
		int last_num = pvt[m].nthr - 1;
		double thresh = pvt[m].thr[last_num];
		for (ori = 0; ori < 2; ori++)
		{
			for (k = 0; k < lenp; k++)
			{
				double sco = 0;
				for (p = 0; p < pwm[m].len; p++)
				{
					switch (str[ori][k + p]) {
					case 'a': {sco += pwm[m].wes[p][0]; break; }
					case 'c': {sco += pwm[m].wes[p][1]; break; }
					case 'g': {sco += pwm[m].wes[p][2]; break; }
					case 't': {sco += pwm[m].wes[p][3]; break; }
					default: {break; }
					}
				}
				sco -= pwm[m].min;
				sco /= pwm[m].raz;
				if (m == m_positive)
				{
					mpos_score[n_pos_cases++] = sco;
				}
				if (sco <= thresh)continue;
				if (sco > best_sco)
				{
					best_sco = sco;
					strncpy(best_seq, &str0[ori][k], pwm[m].len);
					best_seq[pwm[m].len] = '\0';
					if (ori == 0)
					{
						best_ori = '+';
						best_pos = k + 1;
					}
					else
					{
						best_ori = '-';
						best_pos = lenp - k - 1;
					}
				}
			}
		}
		if (best_sco < thresh)
		{
			/*if(m==m_positive)
			{
				a->fit = 0;
				return a->fit;
			}*/
			rfit[m] = 0;
			continue;
		}
		if (m != m_positive)
		{
			if (best_sco >= pvt[m].thr[0])
			{
				rfit[m] = pvt[m].lgpv[0];
				if (ret < rfit[m])
				{
					ret = rfit[m];
					//a->bpos = best_pos;
					a->bpos = best_pos + pwm[m].len / 2;
					a->bmat = m;
					a->ori = best_ori;
					strncpy(a->seq, best_seq, pwm[m].len);
					a->seq[pwm[m].len] = '\0';
					continue;
				}
			}
			else
			{
				for (k = last_num; k >= 1; k--)
				{
					if (best_sco >= pvt[m].thr[k] && best_sco < pvt[m].thr[k - 1])
					{
						rfit[m] = pvt[m].lgpv[k];
						if (ret < rfit[m])
						{
							ret = rfit[m];
							a->bpos = best_pos + pwm[m].len / 2;
							a->bmat = m;
							a->ori = best_ori;
							strncpy(a->seq, best_seq, pwm[m].len);
							a->seq[pwm[m].len] = '\0';
							break;
						}
					}
				}
			}
		}
		else rfit[m] = 0;
	}
	if ((m_positive >= 0 && m_positive < n_pwm) && n_pos_cases < n_pos_sites)return 0;
	a->bpv_n = ret;
	/*	for(m=0;m<n_pwm;m++)
	{
		printf("%.3f\t",rfit[m]);
		if((m+1)%16==0)printf("\n");
	}*/
	qsort(rfit, n_pwm, sizeof(double), compare_qq);	
	/*{
		for (m = 0; m < n_pwm; m++)
		{
			printf("%.3f\t", rfit[m]);
			if (rfit[m] == 0)break;
			if ((m + 1) % 16 == 0)printf("\n");
		}
	}*/
	double wei_n = 0, wei_y = 0;
	for (i = 0; i < 25; i++)wei_n += rfit[i] * wfit[i];
	wei_n /= sum_rfit;
	a->bpv_wn = wei_n;
	if (m_positive < 0 || m_positive >= n_pwm)a->bpv_wy = 1;
	else
	{
		qsort(mpos_score, n_pos_cases, sizeof(double), compare_qq);
		for (i = 0; i < n_pos_sites; i++)wei_y += mpos_score[i];
		wei_y /= n_pos_sites;
		if (wei_y >= pvt[m_positive].thr[0])a->bpv_wy = pvt[m_positive].lgpv[0];
		else
		{
			int last_num = pvt[m_positive].nthr - 1;
			if (wei_y < pvt[m_positive].thr[last_num])a->bpv_wy = pvt[m_positive].lgpv[last_num];
			else
			{
				for (k = last_num; k >= 1; k--)
				{
					if (wei_y >= pvt[m_positive].thr[k] && wei_y < pvt[m_positive].thr[k - 1])
					{
						a->bpv_wy = pvt[m_positive].lgpv[k];
						break;
					}
				}
			}
		}
	}
	//a->fit = a->bpv_wn;
	a->fit = a->bpv_wn / a->bpv_wy;
	//if (a->bpv_wy == 0) {a->fit = 0;}
	//else a->fit = a->bpv_wn / a->bpv_wy;
	//printf("F %f P %f N %f\n", a->fit, a->bpv_wy, a->bpv_wn); 
	return a->fit;
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
int Mut10(town *a, int size, int max, int rda)
{
	int i;
	int rnet, o_da, o_net;
	int sizem=max-size;
	//rda=rand()%size;	
	o_da=a->da[rda];
	rnet = rand()%sizem;
	int ri=0;
	o_net=-1;
	for(i=0;i<max;i++)
	{
		if(a->deg[i]==0)
		{
			if(ri==rnet)
			{
				o_net=i;
				break;
			}
			ri++;
		}
	}
	if(o_net==-1)
	{
		printf("Mutation error!\n");
		return -1;
	}
	int buf=a->da[rda];
	a->da[rda]=o_net;
	a->deg[o_da]--;
	a->deg[o_net]++;
	for(i=0;i<size;i++)
	{
		int d=a->deg[i];
		if(d<0 || d>1)
		{
			printf("Mutation Deg error! %d\n",d);
			exit(1);
		}
	}
	return 1;
}
int Mut11(town *a, int size, int r1)
{
	int r2;	
	//r1=rand()%size;	
	r2=rand()%(size-1);	
	if(r2>=r1)r2++;
	int buf=a->da[r1];
	a->da[r1]=a->da[r2];
	a->da[r2]=buf;
	return 1;
}
int Reco2(town *a, town *b, int size, int max, int &ra, int &rb)
{
	int i,j;
	int x1, x2, y1, y2;	
	int size1=size-1, sizem=max-size-1;
	int inx[NOLIG];
	for(i=0;i<size1;i++)inx[i]=i;		
	ra=rand()%(size1);	
	x1=a->da[ra], x2=a->da[ra+1];
	BigMix1(inx,size1);
	for(i=0;i<size1;i++)
	{
		y1=b->da[inx[i]], y2=b->da[inx[i]+1];
		if(x1==y1 || x2==y2)
		{
			if(x1==y1 && x2==y2){return -1;}
			{rb=inx[i];break;}
		}
	}				
	if(rb==-1)return -1;
	int buf, rb1=rb+1, ra1=ra+1;
	if(x1==y1)
	{
		a->deg[x2]--;
	//	printf("X1Y1 AX2 %d\n",a->deg[x2]);
		a->deg[y2]++;
		b->deg[y2]--;
	//	printf("X1Y1 BY2 %d\n",b->deg[y2]);
		b->deg[x2]++;
		if(a->deg[y2]==2)
		{
			a->deg[y2]=1;			
			for(i=0;i<size;i++)
			{
				if(i==ra || i==ra1)continue;
				if(a->da[i]==y2)
				{
					int r=rand()%sizem;
					if(r==x2)r++;
					int mut=-1, zero=0;
					for(j=0;j<max;j++)
					{
						if(a->deg[j]==0)
						{
							if(zero==r){mut=j;break;}
							zero++;
						}
					}
					if(mut==-1)return -1;
					a->da[i]=mut;
					a->deg[mut]=1;
					break;
				}
			}
		}
		if(b->deg[x2]==2)
		{
			b->deg[x2]=1;			
			for(i=0;i<size;i++)
			{
				if(i==rb || i==rb1)continue;
				if(b->da[i]==x2)
				{
					int r=rand()%sizem;
					if(r==x2)r++;
					int mut=-1, zero=0;
					for(j=0;j<max;j++)
					{
						if(b->deg[j]==0)
						{
							if(zero==r){mut=j;break;}
							zero++;
						}
					}
					if(mut==-1)return -1;
					b->da[i]=mut;
					b->deg[mut]=1;
					break;
				}
			}
		}
		buf=x2;
		a->da[ra1]=y2;//x2=y2;
		b->da[rb1]=buf;//y2=buf;
	}
	else//x2==y2
	{
		a->deg[x1]--;
	//	printf("X2Y2 AX1 %d\n",a->deg[x1]);
		a->deg[y1]++;
		b->deg[y1]--;
	//	printf("X2Y2 BY1 %d\n",b->deg[y1]);
		b->deg[x1]++;
		if(a->deg[y1]==2)
		{
			a->deg[y1]=1;			
			for(i=0;i<size;i++)
			{
				if(i==ra || i==ra1)continue;
				if(a->da[i]==y1)
				{
					int r=rand()%sizem;
					if(r==x2)r++;
					int mut=-1, zero=0;
					for(j=0;j<max;j++)
					{
						if(a->deg[j]==0)
						{
							if(zero==r){mut=j;break;}
							zero++;
						}
					}
					if(mut==-1)return -1;
					a->da[i]=mut;
					a->deg[mut]=1;
					break;
				}
			}
		}
		if(b->deg[x1]==2)
		{
			b->deg[x1]=1;			
			for(i=0;i<size;i++)
			{
				if(i==rb || i==rb1)continue;
				if(b->da[i]==x2)
				{
					int r=rand()%sizem;
					if(r==x2)r++;
					int mut=-1, zero=0;					
					for(j=0;j<max;j++)
					{
						if(b->deg[j]==0)
						{
							if(zero==r){mut=j;break;}
							zero++;
						}
					}
					if(mut==-1)return -1;
					b->da[i]=mut;
					b->deg[mut]=1;
					break;
				}
			}
		}
		buf=x1;
		a->da[ra]=y1;
		b->da[rb]=buf;
	}	
	return 1;
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
	int i, j, k, m;
	char str[1000], filei_seq[300], filei_mat[300], fileo[300], file_log[300], path_pos[300], path_neg[300];
	char mname[50], mext[] = ".pwm", mpv_ext[] = ".dist";//_ups1500_sco_pval
	FILE *in, *out, *outlog;
	if (argc != 9)
	{
		puts("Syntax: 1path pos_matrix 2path neg_matrices 3file input_sequence 4int size(4islo monomerov) 5int matrix_count 6char matrix_name 7file out 8file_log");
		return -1;
	}
	strcpy(path_pos, argv[1]);
	strcpy(path_neg, argv[2]); 
	strcpy(filei_seq, argv[3]);	
	int n_olig = atoi(argv[4]);//4islo monomerov v cepo4ke
	int n_pwm = atoi(argv[5]);
	strcpy(mname, argv[6]);
	strcpy(fileo, argv[7]);
	strcpy(file_log, argv[8]);
	int m_positive = 0;
	int n_pos_sites = 2 * n_olig;//arf
	//int m_ignore[]={ 2, 82, 83, 84, 102, 183, 184, 186, 207, 212, 213, 217, 286, 299, 300, 303, 439, -1 };// 2020 from mcot, 82,83,84 = ARFs
//	int m_ignore[]={ 2, 102, 183, 184, 186, 207, 212, 213, 217, 227, 286, 299, 300, 303, 439, -1 };// 2020 from mcot, 227 = EIN3
//	int m_ignore[] = { 2, 4, 102, 183, 184, 186, 207, 212, 213, 217, 227, 286, 299, 300, 303, 439, -1 };// 2023, 4 = FUS3, 227 = EIN3
	//int m_ignore[] = { 2, 34, 35, 36, 37, 41, 75, 81, 102, 183, 184, 186, 207, 212, 213, 217, 286, 299, 300, 303, 439, -1 };//2021 cbf3
	//int m_ignore[] = { 2, 102, 183, 184, 186, 207, 212, 213, 214, 217, 229, 286, 299, 300, 303, 439, -1 };//2021 x = 214 camta1,229 far
//	int m_ignore[] = { 2, 4, 102, 183, 184, 186, 207, 212, 213, 217, 227, 286, 299, 300, 303, 439, -1 };// 2024, 4 = FUS3, 227 = EIN3
//	int m_ignore[] = { 2, 102, 183, 184, 186, 207, 212, 213, 217, 227, 286, 299, 300, 303, 439, -1 };// 2024, 227 = EIN3
	int m_ignore[] = { 2, 102, 183, 184, 186, 207, 212, 213, 217, 227, 286, 299, 300, 303, 439, -1 };// no anchor 2024

	int m_ignore_count=0;
	for(i=0;m_ignore[i]!=-1;i++)m_ignore_count++;	
//	for(i=0;m_ignore[i]!=-1;i++)m_ignore[i]--;

//oligi
	if((in=fopen(filei_seq,"rt"))==NULL)
	{
		 printf("Input file %s can't be opened!",filei_seq);
		 return -1;
	}
	fgets(str,sizeof(str),in);//head
	char sym=str[0];
	fgets(str,sizeof(str),in);
	DelChar(str,'\n');
	DelChar(str, '\r');
	int olen; //dlina monomera
	int max_olig=0;//4islo monomerov dostupnoe dlya vybora	
	olen=strlen(str);
	rewind(in);
	while(fgets(str,sizeof(str),in)!=NULL)
	{
		if(str[0]==sym)max_olig++;
	}
	rewind(in);
	seq = new char*[max_olig];				
	if(seq==NULL){puts("Out of memory...");exit(1);}
	for(i=0;i<max_olig;i++)
	{			
		seq[i] = new char[olen+1];
		if(seq[i]==NULL){puts("Out of memory...");exit(1);}
	}	
	rewind(in);
	seq0 = new char*[max_olig];				
	if(seq0==NULL){puts("Out of memory...");exit(1);}
	for(i=0;i<max_olig;i++)
	{			
		seq0[i] = new char[olen+1];
		if(seq0[i]==NULL){puts("Out of memory...");exit(1);}
	}	
	i=0;
	while(fgets(str,sizeof(str),in)!=NULL)
	{
		if(strchr("atgcATGC",str[0])!=NULL)
		{
			DelChar(str,'\n');
			DelChar(str, '\r');
			strcpy(seq0[i],str);
			strcpy(seq[i],str);
			TransStr(seq[i]);
			i++;
		}
	}
	fclose(in);
	char sep = '\t';
//matricy
	pwm = new matr[n_pwm];
	if(pwm==NULL ){puts("Out of memory...");return -1;}
	for(m=0;m<n_pwm;m++)
	{
		if(m==0)strcpy(filei_mat, path_pos);
		else strcpy(filei_mat, path_neg);
		strcat(filei_mat,mname);
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
					int u=UnderStol(str,j,sco, sizeof(sco),sep);
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
		strcat(filei_mat,mname);
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
		double pv_thr1 = 5, pv_thr2 = 3;//pv_thr1=1E-5, pv_thr2=0.001;
		while(fgets(str,sizeof(str),in)!=NULL)
		{	
			char sco[50];
			if(isdigit(str[0]))
			{
				int u=UnderStol(str,1,sco,sizeof(sco), sep);
				if(u==-1)
				{				
					printf("Error reading %s file\n",filei_mat);
					return -1;
				}
				double pv=atof(sco);
				//pv = pow(double(10), -pv);
				if(pv>pv_thr1)continue;
				if(pv>pv_thr2)n_str++;
				else 
				{
					n_str++;
					break;
				}
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
		while(fgets(str,sizeof(str),in)!=NULL)
		{	
			char sco[50];
			if(isdigit(str[0]))
			{
				int u=UnderStol(str,1,sco, sizeof(sco), sep);
				if(u==-1)
				{				
					printf("Error reading %s file\n",filei_mat);
					return -1;
				}
				double pv=atof(sco);
				//pv = pow(double(10), -pv);
				if(pv<=pv_thr1)
				{
					pvt[m].lgpv[0] = pv;// - log10(pv);
					pvt[m].thr[0]=atof(str);				
					break;
				}
			}
		}		
		for(i=1;i<pvt[m].nthr;i++)		
		{
			if(fgets(str,sizeof(str),in)==NULL)			
			{
				printf("Unexpected end of %s file\n",filei_mat);
				return -1;
			}
			if(isdigit(str[0]))
			{
				char sco[50];				
				int u=UnderStol(str,1,sco, sizeof(sco), sep);
				if(u==-1)
				{				
					printf("Error reading %s file\n",filei_mat);
					return -1;
				}
				double pv=atof(sco);			
				pvt[m].lgpv[i] = pv;// -log10(pv);
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
	srand( (unsigned)time( NULL ) );
	int m_success_no_max=5, r_success_no_max=1, total_m_success, total_r_success;
	for(i=0;i<POOL;i++)
	{		
	//	printf("\b\b\b%3d",i);
		int gom=0;
		do
		{
			dmut.init(n_olig,max_olig);
			//dmut.print_all(n_mut,n_olig);
			for(j=0;j<i;j++)
			{
				if(GomTown(dmut,pop[j],n_olig)==1)
				{
					gom=1;
					break;
				}
			}
		}
		while(gom==1);
		//dmut.print_all(n_olig,max_olig);
		dmut.get_copy(&pop[i],n_olig,max_olig);		
	//	EvalFit(&pop[i],olen,n_olig,mname,n_pwm,mpv_ext,m_ignore, m_ignore_count);	
		EvalFit(&pop[i], olen,n_olig, mname, n_pwm, mpv_ext, m_positive, n_pos_sites, m_ignore);
	}
	qsort((void*)pop,POOL, sizeof(pop[0]), compare_pop);	
	if ((outlog = fopen(file_log, "wt")) == NULL)
	{
		fprintf(outlog, "Input file %s can't be opened!\n", file_log);
		exit(1);
	}
	for(i=ELIT-1;i>=0;i--)pop[i].print_all(n_olig,max_olig,outlog);	
	for(i=0;i<POOL;i++)
	{
		int check=pop[i].check(n_olig,max_olig);
		if(check==-1)
		{
			printf("After Init\n");
			pop[i].print_all(n_olig,max_olig,outlog);
			exit(1);
		}
	}	
	int iter=0;
	double av_elit_score_prev=0;
	for(i=0;i<ELIT;i++)av_elit_score_prev+=pop[i].fit;
	av_elit_score_prev/=ELIT;
	do
	{
		total_m_success=0;
		total_r_success=0;
		fprintf(outlog,"Mutation %d    ",iter+1);
		for(i=0;i<POOL;i++)
		{
			int m_success_no=0;
		//	if((i+1)%10==0)printf("\b\b\b%3d",i+1);
		//	printf("%3d\n", i + 1);
			if(pop[i].fit==0)continue;
			do
			{
				pop[i].get_copy(&dmut,n_olig,max_olig);
				int mut, mut_type=rand()%2, rda;
				if(dmut.bmat==0)rda=rand()%n_olig;
				else
				{
					int r_pos=int(dmut.bpos/olen);
					int rmin=Max(0,n_olig-1);
					int rmax=Min(n_olig-1, rmin+1);
					int dr=rmax-rmin+1;
					int rdr=rand()%dr;
					rda=rmin+rdr;
				}
				if(mut_type==0)mut=Mut10(&dmut,n_olig,max_olig,rda);		
				else mut=Mut11(&dmut,n_olig,rda);
				if(mut==-1 || dmut.check(n_olig,max_olig)==-1)
				{
					m_success_no++;
					continue;
				}
				int gom=0;
				for(j=0;j<POOL;j++)
				{
					if(GomTown(dmut,pop[j],n_olig)==1){gom=1;break;}
				}
				if(gom==1)
				{
					m_success_no++;
					continue;
				}
				//	dmut.print_all(size,max);				
			//	EvalFit(&dmut,olen,n_olig,mname,n_pwm,mpv_ext,m_ignore,m_ignore_count);				
				EvalFit(&dmut, olen, n_olig, mname, n_pwm, mpv_ext, m_positive, n_pos_sites, m_ignore);
				if(dmut.fit<pop[i].fit)
				{
					dmut.get_copy(&pop[i],n_olig,max_olig);	
					m_success_no=0;
					total_m_success++;
					int check=pop[i].check(n_olig,max_olig);
					if(check==-1)
					{
						fprintf(outlog,"After Mut!\n");
						exit(1);
					}
				}
				else m_success_no++;
			}
			while(m_success_no<m_success_no_max);
		}
		qsort((void*)pop,POOL, sizeof(pop[0]), compare_pop);	
		fprintf(outlog,"\nGA %d\tMut %d\n",iter+1,total_m_success);	
		for(i=ELIT-1;i>=0;i--)pop[i].print_all(n_olig,max_olig,outlog);
		fprintf(outlog,"Recombination %d      ",iter+1);
		int olen32=olen+olen/2;
		int pool_last=POOL-1;
		for(i=pool_last;i>0;i--)
		{			
			//if((i+1)%10==0)			
			if(pop[i].fit==0)continue;
		//	printf("Rec 1st %3d\n",i);			
			for(j=i-1;j>=0;j--)
			{
				//if((j+1)%50==0)printf("\b\b\b\b\b\b\b%3d %3d",i,j);
			//	if(i==j)continue;
				if(pop[j].fit==0)continue;
			//	printf("Rec 2nd %3d\n",j);						
				int r_success_no=0;
				do
				{
					int rec;
					pop[i].get_copy(&drec[0],n_olig,max_olig);			
					pop[j].get_copy(&drec[1],n_olig,max_olig);	
					int ra=-1, rb=-1;
			//		printf("\n");
			//		for(k=0;k<2;k++)drec[k].print_all(n_olig,max_olig);									
					if(drec[0].bpos<=olen32)ra=0;
					else
					{
						ra=1+(drec[0].bpos-olen32)/olen;
					}
					//for(k=0;k<2;k++)drec[k].print_all(n_olig,max_olig);
					rec = Reco2(&drec[0],&drec[1],n_olig,max_olig,ra,rb);
			//		for(k=0;k<2;k++)drec[k].print_all(n_olig,max_olig);
					if((rec==-1 || (drec[0].check(n_olig,max_olig)==-1 || drec[1].check(n_olig,max_olig)==-1)))					
					{
						r_success_no++;
						continue;
					}					
					int gom=0;
					for(k=0;k<2;k++)
					{						
						for(m=0;m<POOL;m++)
						{
							if(GomTown(drec[k],pop[m],n_olig)==1){gom=1;break;}
						}
						if(gom==1)break;
					}
					if(gom==1)
					{
						r_success_no++;
						break;
					}
					double fitmin=Min(pop[i].fit,pop[j].fit);
					for(k=0;k<2;k++)EvalFit(&drec[k], olen, n_olig, mname, n_pwm, mpv_ext, m_positive, n_pos_sites, m_ignore);
					if(drec[0].fit<fitmin || drec[1].fit<fitmin)
					{
					//	for(k=0;k<2;k++)drec[k].print_all(size,max);	
						drec[0].get_copy(&pop[i],n_olig,max_olig);			
						drec[1].get_copy(&pop[j],n_olig,max_olig);						
						total_r_success++;
					}
					else r_success_no++;										
				}
				while(r_success_no<r_success_no_max);
			}
		}				
		qsort((void*)pop,POOL, sizeof(pop[0]), compare_pop);
		for(i=ELIT-1;i>=0;i--)pop[i].print_all(n_olig,max_olig,outlog);
		iter++;
		double av_elit_score=0;
		for(i=0;i<ELIT;i++)av_elit_score+=pop[i].fit;
		av_elit_score/=ELIT;
		fprintf(outlog,"\nGA %d\tMut %d\tRec %d\tCurFitElit %f\tPrevFitElit %f\n",iter+1,total_m_success,total_r_success,av_elit_score,av_elit_score_prev);	
		fclose(outlog);
		if ((outlog = fopen(file_log, "at")) == NULL)
		{
			fprintf(outlog, "Input file %s can't be opened!\n", file_log);
			exit(1);
		}
		if(av_elit_score>=av_elit_score_prev)break;
		av_elit_score_prev=av_elit_score;		
	}
	while(total_m_success+total_r_success>POOL);
	fprintf(outlog,"Final\n");
	fclose(outlog);
	//int ext_elit= Min(POOL, 5*ELIT);
	//for(i=0;i<ext_elit;i++)pop[i].print_all(n_olig,max_olig);	
	if((out=fopen(fileo,"wt"))==NULL)
	{
		printf("Output file %s can't be opened!\n", fileo);
		return -1;
	}	
	for(k=0;k<POOL;k++)
	{
		fprintf(out, ">%d\t", k + 1);
		fprintf(out, "Fit %.6f\tPvNeg %.6f\tPvNegW %.6f\tPvPosW %.6f\t", pop[k].fit, pop[k].bpv_n, pop[k].bpv_wn, pop[k].bpv_wy);
		fprintf(out, "M %d\tP %d%c %s", pop[k].bmat, pop[k].bpos, pop[k].ori, pop[k].seq);
		for(i=0;i<n_olig;i++)fprintf(out,"\t%d",1+pop[k].da[i]);
		fprintf(out, "\n");
		for(i=0;i<n_olig;i++)fprintf(out,"%s",seq0[pop[k].da[i]]);
		fprintf(out,"\n");
	}
	fclose(out);
}
