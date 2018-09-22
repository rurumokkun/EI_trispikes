///// 3つのグループ間のインタラクション   /////

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
using namespace std;
#include <stdlib.h> 
#include <time.h>
//#include <cstdlib>

////  '12ichi
//  new	to Wspikenohikaku
//  12/8 グループ2の重みのヒストグラム作成と周波数解析を行う予定
#define getrandom(max1) ((rand()%(int)((max1)))) // random integer between 0 and max-1     //getrandom(N)は0~N-1のいずれかを出力
#define NM 256
const   int		Ngr=3;                  	//ニューロングループの数 
const   int		SELF_SPIKE_TIME =1000;	                                //グループ別自己発火時間
const	int		INTERACT_TIME	=2000;					//グループ間相互作用時間
const	int		SIM_TIME=SELF_SPIKE_TIME+INTERACT_TIME;                //シミュレーション時間(シミュレーション前に決定)
const	int		NUM_EI_PATTERN =9;

const   int             num_trial=20;             
//const	int		Ne = 800;		// excitatory neurons			
//const	int		Ni = 200;		// inhibitory neurons				 
const	int		N  = 1000;		// total number of neurons	
const   int     	interM=3;		// 1ニューロンあたりのnyu-ronnグループ間シナプス数		
const	int		M  = 100;		// the number of synapses per neuron 　
const	int		D  = 20;		// maximal axonal conduction delays
//const	int		maxE=810;

float	sm = 10.0;		// maximal synaptic strength		
int		post[Ngr][N][M];				// indeces of postsynaptic neurons
float	s[Ngr][N][M], sd[Ngr][N][M];		// matrix of synaptic weights and their derivatives
short	delays_length[Ngr][N][D];	// distribution of delays
short	delays[Ngr][N][D][M];		// arrangement of delays   
int		N_pre[Ngr][N], I_pre[Ngr][N][3*M], D_pre[Ngr][N][3*M];	// presynaptic information
int		interN_pre[Ngr][Ngr-1][N], interI_pre[Ngr][Ngr-1][N][10*interM], interD_pre[Ngr][Ngr-1][N][10*interM];
float	*s_pre[Ngr][N][3*M], *sd_pre[Ngr][N][3*M];		// presynaptic weights
float	*inter_s_pre[Ngr][Ngr-1][N][30*interM], *inter_sd_pre[Ngr][Ngr-1][N][30*interM];		// グループ間シナプス前結合荷重&更新用変数
float	LTP[Ngr][N][1001+D], LTD[Ngr][N];	// STDP functions : LTP(Long Term Potentiation) LTD(Long Term Potentiation)
float	inter_LTP[Ngr][N][1001+D+10],inter_LTD[Ngr][N];
float	a[Ngr][N], d[Ngr][N];				// neuronal dynamics parameters
float	v[Ngr][N], u[Ngr][N];				// activity variables
int		N_firings[Ngr];				// the number of fired neurons 
const int N_firings_max=100*N*10;	// upper limit on the number of fired neurons per sec		//発火数の上限: 変えていいのか？論文を確認すべき 9:1は必ず超える
int		firings[Ngr][N_firings_max][2]; // indeces and timings of spikes
int         before_firings[Ngr][N_firings_max][2];  	// 
int         interpost[Ngr][Ngr-1][N][interM];		//
float       interpost_s[Ngr][Ngr-1][N][interM];	//グループ間結合重み
float       interpost_sd[Ngr][Ngr-1][N][interM];	//グループ間結合重みの更新用変数
float       interpost_delays[Ngr][Ngr-1][N][interM];	//グループ間伝達遅れ
float       I[Ngr][N];				// 入力
int		z[Ngr];				// N_firingsコピー用
int		y[Ngr];

int	EI_1,EI_2,EI_3;
int num;  //E/I組み合わせ指定用変数

int               Ne[Ngr][NUM_EI_PATTERN]={{850,840,830,820,810,800,700,600,500},     {850,840,830,820,810,800,700,600,500},     {850,840,830,820,810,800,700,600,500}};

int                 Ni[Ngr][NUM_EI_PATTERN];   //1/20実験 真ん中変える

///////////////////////////
///	grandagrtoogr   ///   //グループとそれの指し示す方向の写像=>グループ
///////////////////////////
int grandagrtoogr(int gr, int agr)
{
     int ogr;

   if(gr==0&&agr==0)
	ogr=1;
 else if(gr==1&&agr==1)
	ogr=0;
 else if(gr==1&&agr==0)
 	ogr=2;
 else if(gr==2&&agr==1)
  	ogr=1;
 else if(gr==2&&agr==0)
   	ogr=0;
 else   ogr=2;


    return(ogr);
}


///////////////////////////
//    agrcahnge		///
///////////////////////////
int agrchange(int agr)
{

   int ogr;
 
   ogr=abs(agr-1);

   return(ogr);
　
}
///////////////////////////
///       初期化	///
///////////////////////////
void initialize()// initialize　start				
｛	int i,j,k,jj,dd, exists,gr, r ; 		//gr:ニューロングループ指定用変数
        
       
	for(gr=0;gr<Ngr;gr++){　　//for1

         for (i=0;i<Ne[gr][num];i++) a[gr][i]=0.02;// RS type
     
	for (i=Ne[gr][num];i<N;i++) a[gr][i]=0.1;  // FS type
    
	for (i=0;i<Ne[gr][num];i++) d[gr][i]=8.0;  // RS type
	for (i=Ne[gr][num];i<N;i++) d[gr][i]=2.0;  // FS type
	srand((unsigned) time(NULL));
	for (i=0;i<N;i++) for (j=0;j<M;j++) 
	{//for2
		do{//do1
			exists = 0;		// avoid multiple synapses
			if (i<Ne[gr][num]) r = getrandom(N);
			else	  r = getrandom(Ne[gr][num]);// inh -> exc only
			if (r==i) exists=1;									// no self-synapses 
			for (k=0;k<j;k++) if (post[gr][i][k]==r) exists = 1;	// synapse already exists  
		}while (exists == 1);//do1　end
		post[gr][i][j]=r;
	}//for2　end
	for (i=0;i<Ne[gr][num];i++)	for (j=0;j<M;j++) s[gr][i][j]=6.0;  // initial exc. 興奮性シナプス荷重
	for (i=Ne[gr][num];i<N;i++)	for (j=0;j<M;j++) s[gr][i][j]=-5.0; // 		    抑制性シナプス荷重
  	for (i=0;i<N;i++)	for (j=0;j<M;j++) sd[gr][i][j]=0.0; // synaptic derivatives 
  	for (i=0;i<N;i++) 	     		//delays, delays_lengthの初期化:  興奮性,抑制性で場合分け
	{//for3
		short ind=0;
		if (i<Ne[gr][num])
		{//if　1
			for (j=0;j<D;j++) //for4
			{	delays_length[gr][i][j]=5;	//(M/D=5) uniform distribution of exc. synaptic delays
				for (k=0;k<delays_length[gr][i][j];k++)  // 0<=k<5
					delays[gr][i][j][k]=ind++;
			}//for4　end
		}//if　1　end
		else
		{//else　1
			for (j=0;j<D;j++) delays_length[gr][i][j]=0;
			delays_length[gr][i][0]=M;			// all inhibitory delays are 1 ms
			for (k=0;k<delays_length[gr][i][0];k++)
					delays[gr][i][0][k]=ind++;
		}//else1　end
	}//　for3　end
	
  	for (i=0;i<N;i++)
	{//for5
		N_pre[gr][i]=0;
		for (j=0;j<Ne[gr][num];j++)
		for (k=0;k<M;k++)
		if (post[gr][j][k] == i)		// find all presynaptic neurons 
		{
			I_pre[gr][i][N_pre[gr][i]]=j;	// add this neuron to the list
			for (dd=0;dd<D;dd++)	// find the delay
				for (jj=0;jj<delays_length[gr][j][dd];jj++)
					if (post[gr][j][delays[gr][j][dd][jj]]==i) D_pre[gr][i][N_pre[gr][i]]=dd;
			s_pre[gr][i][N_pre[gr][i]]=&s[gr][j][k];	// pointer to the synaptic weight	
			sd_pre[gr][i][N_pre[gr][i]++]=&sd[gr][j][k];// pointer to the derivative
		}
	}//for5　end

	for (i=0;i<N;i++)	for (j=0;j<1+D;j++) LTP[gr][i][j]=0.0;
	for (i=0;i<N;i++)	LTD[gr][i]=0.0;
	for (i=0;i<N;i++)	v[gr][i]=-65.0;		// initial values for v
	for (i=0;i<N;i++)	u[gr][i]=0.2*v[gr][i];	// initial values for u

	N_firings[gr]=1;		// spike timings
	firings[gr][0][0]=-D;	// put a dummy spike at -D for simulation efficiency 
	firings[gr][0][1]=0;	// index of the dummy spike
    
  	before_firings[gr][0][0]=-D;  
        before_firings[gr][0][1]=0;																																																																																																																																																																																																																																																																																														
    }	//grについてのfor1　end
  }   //initialize関数の終わり





///////////////////////////
///	グループ初期化	///
///////////////////////////
void group_initialize()                 
{//　group　initialize　start
  int gr,agr,i,j,r,k,exists;   //agr:=another group


   for(gr=0;gr<Ngr;gr++){//for6　
    for(agr=0;agr<Ngr-1;agr++){//for7
    for(i=0;i<N;i++){//for8
     for(j=0;j<interM;j++){//for9
      				
       do{
         exists=0;
        if(i<Ne[gr][num])r=getrandom(N);				
        else    r=getrandom(Ne[gr][num]);
    	for(k=0;k<j;k++)if (interpost[gr][agr][i][k]==r)exists=1;		
   
        }while(exists==1);



        interpost[gr][agr][i][j]=r;					//シナプス後ニューロン
        if(i<Ne[gr][num] && (gr==0&&agr==0 ||gr==1||gr==2&&agr==1))  interpost_s[gr][agr][i][j]=6.0;			
        else if(i>=Ne[gr][num] && (gr==0&&agr==0 ||gr==1||gr==2&&agr==1))    interpost_s[gr][agr][i][j]=-5.0;
	else interpost_s[gr][agr][i][j]=0;	
	interpost_sd[gr][agr][i][j]=0;					//11/08追加
   	interpost_delays[gr][agr][i][j]=getrandom(21)+10.0;			//伝達遅れの初期化: 10~30の間
      }//jに関するforループ　//for9　end
    }//iに関するforループ　//for8　end
 }//for7　end
}//for6　end

//////////////		グループ間の種々の変数の決定		/////////////////////////////////////////////////////////////
	for(gr=0;gr<Ngr;gr++)
        for(agr=0;agr<Ngr-1;agr++)
	for(i=0;i<N;i++)//for10　
	{
		interN_pre[gr][agr][i]=0;
		for(j=0;j<Ne[gr][num];j++)
		for(k=0;k<interM;k++)
		if(interpost[gr][agr][j][k]==i)
		{
			interI_pre[grandagrtoogr(gr,agr)][agrchange(agr)][i][interN_pre[grandagrtoogr(gr,agr)][agrchange(agr)][i]]=j;
		
		
			interD_pre[grandagrtoogr(gr,agr)][agrchange(agr)][i][interN_pre[grandagrtoogr(gr,agr)][agrchange(agr)][i]]=interpost_delays[gr][agr][j][k];

			inter_s_pre[grandagrtoogr(gr,agr)][agrchange(agr)][i][interN_pre[grandagrtoogr(gr,agr)][agrchange(agr)][i]]=&interpost_s[gr][agr][j][k];
			inter_sd_pre[grandagrtoogr(gr,agr)][agrchange(agr)][i][interN_pre[grandagrtoogr(gr,agr)][agrchange(agr)][i]++]=&interpost_sd[gr][agr][j][k];
		}

	}//for10　end
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  }
//}


/////////////   2/1追加	/////////////////////////////////////////////////////////////////////////////////////////////
	for(gr=0;gr<Ngr;gr++)for (i=0;i<N;i++)	for (j=0;j<1+10+D;j++) inter_LTP[gr][i][j]=0.0;   //伝達遅れが10~30それぞれに対して0で初期化
	for(gr=0;gr<Ngr;gr++)for (i=0;i<N;i++)	inter_LTD[gr][i]=0.0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



}//group　initialize終わり




int main()
{      
	
        char fname_spikes_1[50],fname_spikes_2[50],fname_spikes_3[50];
        char fname_LAP_1[50],fname_LAP_2[50],fname_LAP_3[50];
	char fname_interpost_s_1[50],fname_interpost_s_2[50];
	char fname_interpost_s_3[50],fname_interpost_s_4[50];
	char fname_interpost_s_5[50],fname_interpost_s_6[50];
	char fname_s_1[50];
	char fname_s_2[50];
	char fname_s_3[50];
	char str[30];
        char st[N]={'\0'};

	int		i, j, k, sec, t,gr,agr,trial;
       
	int		temp_gr;
        int             l;
        int		ogr;
        int  		temp[Ngr];			//
	time_t		t1,t2,timer;
        struct 		tm *timeptr;
	
	//int 		interpost_s_1_range[101],interpost_s_2_range[101];          //グループ間結合重みのカウント変数:interpost_s_range[k]は重みが0~1の数
	//int 		interpost_s_3_range[101],interpost_s_4_range[101];
	//int 		interpost_s_5_range[101],interpost_s_6_range[101];
	//int		s_range[Ngr][101];			//グループ内結合重みのカウント変数:	s_range[k]
        float		temp_I[Ngr][N];	
	float		v_sum[Ngr],LAP[Ngr][1000];
	
//        float		MI_1,MI_2,MI_3;				//相互情報量
        float           full_LAP[Ngr][1000];//full_LAPはLAPデータの全時間分格納		//secごとにリセット
//	float		LAPx[1000],LAPy[1000],LAPz[1000];
	int            firing_count[NUM_EI_PATTERN][Ngr][N];
        float         firing_rate_seq[NUM_EI_PATTERN][Ngr][SIM_TIME];  
	//float		interpost_s_1_pd[101],interpost_s_2_pd[101];
	//float		interpost_s_3_pd[101],interpost_s_4_pd[101];
	//float		interpost_s_5_pd[101],interpost_s_6_pd[101];
	//float		          s_pd[Ngr][101];
	char	flag;
	int	temp_Ne[Ngr][num];				//一時的にNe[gr]を逃がす
        
	FILE	*fs_1,*fs_2,*fs_3,*fs_LAP_1,*fs_LAP_2,*fs_LAP_3;			      //fs_1はグループ１の発火情報, fs_2はグループ2の発火情報を出力する用
	
	FILE	*fp_interpost_s_1,*fp_interpost_s_2,*fp_interpost_s_3,*fp_interpost_s_4,*fp_interpost_s_5,*fp_interpost_s_6;			//グループ間結合重み
	FILE    *fp_s_1,*fp_s_2,*fp_s_3;					//グループ内結合重み
//	FILE    *fp_MI;
 	FILE    *fp_firing_count;
	FILE    *fp_firing_rate;    //17'11/7
	//FILE	*fp_numofNe;				//異常な発火を示すニューロンのインデックスを格納したい


/////////////////////////////////////////////////////////////////////////////////////////// 	

	
	for(num=0;num<NUM_EI_PATTERN;num++)
	{//for11
	  Ni[0][num]=1000-Ne[0][num];
	  Ni[1][num]=1000-Ne[1][num];
	  Ni[2][num]=1000-Ne[2][num];

	}//for11　end

//////////////////////////////////////////////////////////////////////////////////////////


     
	
         

        //use_getrusage();
       
     t1=time(NULL);		//開始時間

for(trial=0;trial<num_trial;trial++)
{//trial文　//for12　
      timer=time(NULL);
         timeptr=localtime(&timer);
         strftime(st,NM,"data%Y%m%d%H%M",timeptr);

         system("mkdir data_s");
	system("mkdir data_interpost_s");
	system("mkdir data_LAP");
	system("mkdir data_spikes");
        system("mkdir data_firing_count");
 
	 system("mkdir data`date '+%Y%m%d%H%M'`");
        
         sprintf(str,"mv data_* %s",st);


      for(num=0;num<NUM_EI_PATTERN;num++)
      {//for13

	EI_1=Ne[0][num]*10+Ni[0][num]/10;
	EI_2=Ne[1][num]*10+Ni[1][num]/10;
	EI_3=Ne[2][num]*10+Ni[2][num]/10;

        //E_I_initialize(Ne, Ni, k1, k2);
 
        printf("EIバランスは,グループ1が%d:%d\tグループ2が%d:%d\tグループ3が%d:%dです\n",Ne[0][num]/100,Ni[0][num]/100,Ne[1][num]/100,Ni[1][num]/100,Ne[2][num]/100,Ni[2][num]/100);
         	
	initialize();
				// 関数呼び出し: assign connections, weights, etc. 
       group_initialize();	//			　関数呼び出し: グループ初期化
        
	
        for(gr=0;gr<Ngr;gr++)
         for(i=0;i<N;i++)
         firing_count[num][gr][i]=0;                        //発火数を格納
      


       for (sec=0; sec<SIM_TIME; sec++)	// シミュレーション時間ははじめに"SIM_TIME"で指定	
	{//for14	
		for (t=0;t<1000;t++)				// simulation of 1 sec
		{//for15
			for(gr=0;gr<Ngr;gr++){//for16
                          for (i=0;i<N;i++) I[gr][i] = 0.0;	// reset the input 

			if(sec<SIM_TIME/*-100*/){
			for (k=0;k<N/1000;k++)			// k? (k=0のみ?)
				I[gr][getrandom(N)]=20.0;		// random thalamic :input 20の入力   ※ただし,グループ間相互作用時は行わない
			}//for16　end

			for (i=0;i<N;i++) 
			if (v[gr][i]>=30)					// did it fire?
			{
				v[gr][i] = -65.0;					// voltage reset
				u[gr][i]+=d[gr][i];					// recovery variable reset
				LTP[gr][i][t+D]= 0.1;		
				LTD[gr][i]=0.12;
				inter_LTP[gr][i][t+D+10]=0.1/**0.6*/;		//0.1*exp(-10/20)が初期値
				inter_LTD[gr][i]=0.12/**0.6*/;		//0.12*exp(-10/20)が初期値

                                                             
				for (j=0;j<N_pre[gr][i];j++) *sd_pre[gr][i][j]+=LTP[gr][I_pre[gr][i][j]][t+D-D_pre[gr][i][j]-1];// this spike was after pre-synaptic spikes    //sdの更新(pre->postの順に発火)

///////////////////////////////////  inter_sdの更新  /////////////////////////////////////////////////////////////////////////////
				if(sec>=SELF_SPIKE_TIME)
				{
				 for(agr=0;agr<Ngr-1;agr++)
				  //if(gr==0&&agr==0||gr==1||gr==Ngr-1&&agr==1)
                                  {//for17
					for(j=0;j<interN_pre[gr][agr][i];j++) *inter_sd_pre[gr][agr][i][j]+=inter_LTP[gr][interI_pre[gr][agr][i][j]][t+D+10-interD_pre[gr][agr][i][j]];
				   }//for17　end
				   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				firings[gr][N_firings[gr]  ][0]=t;			//[0]: 時間を格納
				firings[gr][N_firings[gr]++][1]=i;			//[1]: ニューロン番号を格納
				if (N_firings[gr] == N_firings_max) {cout << "Two many spikes at t=" << t << " (ignoring all)";N_firings[gr]=1;}
			}
			z[gr]=N_firings[gr];
  			temp[gr]=z[gr];
			while (t-firings[gr][--z[gr]][0] <D)
			{
				for (j=0; j< delays_length[gr][firings[gr][z[gr]][1]][t-firings[gr][z[gr]][0]]; j++)
				{//for18
					i=post[gr][firings[gr][z[gr]][1]][delays[gr][firings[gr][z[gr]][1]][t-firings[gr][z[gr]][0]][j]]; 
					I[gr][i]+=s[gr][firings[gr][z[gr]][1]][delays[gr][firings[gr][z[gr]][1]][t-firings[gr][z[gr]][0]][j]];
					if (firings[gr][z[gr]][1] <Ne[gr][num]) // this spike is before postsynaptic spikes
						sd[gr][firings[gr][z[gr]][1]][delays[gr][firings[gr][z[gr]][1]][t-firings[gr][z[gr]][0]][j]]-=LTD[gr][i];   //sdの更新(post->preの順に発火)
					
					
				}//for18　end
			}
    			


			
                            temp_gr=gr;

///////////////////////////////////////////////	interpost_sdの更新  	////////////////////////////////////////////////////////
	//グループ内のsd更新の後、グループ間のsd更新		///	
			if(sec>=SELF_SPIKE_TIME && gr==Ngr-1)
		     for(gr=0;gr<Ngr;gr++)z[gr]=temp[gr];
	
			for(gr=0;gr<Ngr;gr++)
			for(agr=0;agr<Ngr-1;agr++)
			{ //for19
				
				//if(gr==0&&agr==0||gr==1||gr==Ngr-1&&agr==1)     //グループ1,3間は未接続
  				while(t-firings[grandagrtoogr(gr,agr)][--z[grandagrtoogr(gr,agr)]][0] <D+10 && z[grandagrtoogr(gr,agr)]!=0)
				{
					if(t-firings[grandagrtoogr(gr,agr)][z[grandagrtoogr(gr,agr)]][0]>=10&&firings[grandagrtoogr(gr,agr)][z[grandagrtoogr(gr,agr)]][1]<Ne[grandagrtoogr(gr,agr)][num])

					for(j=0; j< interM; j++)
    					{//for20
					     if(t-firings[grandagrtoogr(gr,agr)][z[grandagrtoogr(gr,agr)]][0]==interpost_delays[gr][agr][firings[grandagrtoogr(gr,agr)][z[grandagrtoogr(gr,agr)]][1]][j])
                                      		{i=interpost[gr][agr][firings[grandagrtoogr(gr,agr)][z[grandagrtoogr(gr,agr)]][1]][j];
				      		 I[gr][i]+=interpost_s[gr][agr][firings[grandagrtoogr(gr,agr)][z[grandagrtoogr(gr,agr)]][1]][j];
				      		 
					    	 interpost_sd[gr][agr][firings[grandagrtoogr(gr,agr)][z[grandagrtoogr(gr,agr)]][1]][j]-=inter_LTD[gr][i];
						 }
					}//for20　end
				}
			}//for19　end		


			for(gr=0;gr<Ngr;gr++) z[gr]=temp[gr];	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







				gr=temp_gr;
				v_sum[gr]=0;
				
				
				
					
				

				   for(i=0;i<Ne[gr][num];i++){
                        	    v_sum[gr]+=v[gr][i];
				   } 
				
   				
				LAP[gr][t]=v_sum[gr]/float(Ne[gr][num]);	
			        full_LAP[gr][t]=LAP[gr][t];

			


		      /*if(sec>=SELF_SPIKE_TIME){							//グループ内発火終了後
     			//temp_I[Ngr-1][N-1]=I[Ngr-1][N-1];
                        intergroup_input(t,gr,sec);						//グループ間の入力を計算
                        //printf("変化前:%f  変化後:%f\n",temp_I[Ngr-1][N-1],I[Ngr-1][N-1]);  
		       }*/



			if(gr==Ngr-1)for(gr=0;gr<Ngr;gr++)

			for (i=0;i<N;i++)
			{//for21
				v[gr][i]+=0.5*((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]); // for numerical stability
				v[gr][i]+=0.5*((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]); // time step is 0.5 ms
				u[gr][i]+=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
				LTP[gr][i][t+D+1]=0.95*LTP[gr][i][t+D];
				LTD[gr][i]*=0.95;
			}//for21　end


 			gr=temp_gr;	
				//Ne[gr][num]=temp_Ne[gr][num];
			
			/////		グループ間LTP,LTDの更新       //////
                        if(sec>=SELF_SPIKE_TIME &&gr==Ngr-1)for(gr=0;gr<Ngr;gr++)for (i=0;i<N;i++)
			{
				
				inter_LTP[gr][i][t+D+1+10]=0.95*inter_LTP[gr][i][t+D];			//0.95=exp(-1/20)
				inter_LTD[gr][i]*=0.95;
			}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////









                      }//grについてのfor文終わり　//for16　end




		}//tについてのfor文終わり　//for15　end
                 for(gr=0;gr<Ngr;gr++){//for22
                    y[gr]=z[gr];
                   for(i=0;i<z[gr];i++){//for23
                       
          		before_firings[gr][i][0]=firings[gr][i][0];  //time
			before_firings[gr][i][1]=firings[gr][i][1];  //ニューロン番号   
         		}//for23　end
                
                  }//for22　end
 
                for(gr=0;gr<Ngr;gr++){//for24
                 
                /// グループ1ならfs_1に,2ならfs_2に出力  ///
		cout << "sec=" << sec << ", firing rate=" << float(N_firings[gr])/N << ", グループ1=\t" << Ne[0][num]/100<<":"<<Ni[0][num]/100  <<",グループ2=\t" << Ne[1][num]/100<<":"<<Ni[1][num]/100<<",グループ3=\t" << Ne[2][num]/100<<":"<<Ni[2][num]/100<< "\t" <<trial+1<<"回目"<<"\n";	//発火数を出力

　　　　　　　　　　　　　　　　　firing_rate_seq[num][gr][sec]=float(N_firings[gr])/N;

                 if(gr==0){
                 sprintf(fname_spikes_1,"spikes_1_%dsec_%d_%d_%d.txt",SIM_TIME,EI_1,EI_2,EI_3);



                 fs_1 = fopen(fname_spikes_1,"w");
                 

		for (i=1;i<N_firings[gr];i++)
			if (firings[gr][i][0] >=0)
				fprintf(fs_1, "%d  %d\n", firings[gr][i][0], firings[gr][i][1]);
		fclose(fs_1);
                }



                else if(gr ==1){
                 sprintf(fname_spikes_2,"spikes_2_%dsec_%d_%d_%d.txt",SIM_TIME,EI_1,EI_2,EI_3);
                 fs_2 = fopen(fname_spikes_2,"w");


                 for (i=1;i<N_firings[gr];i++)
			if (firings[gr][i][0] >=0)
				fprintf(fs_2, "%d  %d\n", firings[gr][i][0], firings[gr][i][1]);
		fclose(fs_2);

                }


		 else /*gr==2*/{
                 sprintf(fname_spikes_3,"spikes_3_%dsec_%d_%d_%d.txt",SIM_TIME,EI_1,EI_2,EI_3);
                 fs_3 = fopen(fname_spikes_3,"w");


                 for (i=1;i<N_firings[gr];i++)
			if (firings[gr][i][0] >=0)
				fprintf(fs_3, "%d  %d\n", firings[gr][i][0], firings[gr][i][1]);
		fclose(fs_3);

                }





		for (i=0;i<N;i++)		// prepare for the next sec
		 for(agr=0;agr<Ngr-1;agr++)     //'17,1.5追記
			for (j=0;j<D+1;j++)
			LTP[gr][i][j]=LTP[gr][i][1000+j];
		z[gr]=N_firings[gr]-1;

		for (i=0;i<N;i++)		// prepare for the next sec
			for (j=0;j<D+10+1;j++)
			inter_LTP[gr][i][j]=inter_LTP[gr][i][1000+j];

		while (1000-firings[gr][z[gr]][0]<D) z[gr]--;
		for (i=1;i<N_firings[gr];i++)
		{//for25
			firings[gr][i][0]=firings[gr][z[gr]+i][0]-1000;
			firings[gr][i][1]=firings[gr][z[gr]+i][1];
		}//for25　end
		N_firings[gr] = N_firings[gr]-z[gr];

		if(sec<SIM_TIME-200){
		for (i=0;i<Ne[gr][num];i++)	// modify only exc connections:  重みの更新(興奮性シナプスに限る)
		for(agr=0;agr<Ngr-1;agr++)
		for (j=0;j<M;j++)
		{//for26
			s[gr][i][j]+=0.01+sd[gr][i][j];
			sd[gr][i][j]*=0.9;			
			if (s[gr][i][j]>sm) s[gr][i][j]=sm;		//重みが上限10を超えない
			if (s[gr][i][j]<0) s[gr][i][j]=0.0;		//重みが下限0を超えない
		      }//for26　end
	        }                                    

		///グループ間結合重み更新///
		if(sec>=SELF_SPIKE_TIME && sec<SIM_TIME-100 ){
		for (i=0;i<Ne[gr][num];i++)
		for(agr=0;agr<Ngr-1;agr++)
	 	if(gr==0&&agr==0||gr==1||gr==2&&agr==1)
		for (j=0;j<interM;j++)
		{//for27
			interpost_s[gr][agr][i][j]+=0.01+interpost_sd[gr][agr][i][j];
			interpost_sd[gr][agr][i][j]*=0.9;			
			if (interpost_s[gr][agr][i][j]>sm) interpost_s[gr][agr][i][j]=sm;		//重みが上限10を超えない
			if (interpost_s[gr][agr][i][j]<0) interpost_s[gr][agr][i][j]=0.0;		//重みが下限0を超えない
		}//for27　end
			                }
		//////////////////////////

                 }//grについてのfor文終わり　//for24　end




///////////////////////////////
///     LAPの出力        //////
///////////////////////////////
      sprintf(fname_LAP_1,"LAP_1_%d_%d_%d_%dsec.txt",EI_1,EI_2,EI_3,SIM_TIME);
      sprintf(fname_LAP_2,"LAP_2_%d_%d_%d_%dsec.txt",EI_1,EI_2,EI_3,SIM_TIME);
      sprintf(fname_LAP_3,"LAP_3_%d_%d_%d_%dsec.txt",EI_1,EI_2,EI_3,SIM_TIME);
        
      fs_LAP_1=fopen(fname_LAP_1,"a");                                              //"a"は追加書き込み: ファイルがなければ新規作成,あれば末端に追記する
      fs_LAP_2=fopen(fname_LAP_2,"a");
      fs_LAP_3=fopen(fname_LAP_3,"a");

                
                 for(t=0;t<1000;t++) {   //for28                              
                        fprintf(fs_LAP_1,"%f\n",  full_LAP[0][t]);
			fprintf(fs_LAP_2,"%f\n",  full_LAP[1][t]);
			fprintf(fs_LAP_3,"%f\n",  full_LAP[2][t]);
			
                   }//for28　end
	  
                 fclose(fs_LAP_1);
                 fclose(fs_LAP_2);
		 fclose(fs_LAP_3);
                
/////////////////////////////////////////////////////////////////////////////////////////   

		//	if(sec>=SIM_TIME-1){
		//		for(i=0;i<1000;i++)
   		//		{  
                  //         		LAPx[(sec-1999)*1000+i]=full_LAP[0][i];
			//	}

			//	for(i=0;i<1000;i++)
			//	{	
		          // 		LAPy[(sec-1999)*1000+i]=full_LAP[1][i];
				
                        //	}

			//	for(i=0;i<1000;i++)
			//	{	
		          // 		LAPz[(sec-1999)*1000+i]=full_LAP[2][i];
				
                        	//}
			//}
/////////////////////////////////////////////////////////////////////////////////////////

                if(sec==0||sec==1000||sec==1999||sec==2999)
		{
			sprintf(fname_interpost_s_1,"interpost_s_%d_%dsec_%d_%d_%d.txt",1,sec,EI_1,EI_2,EI_3);
			sprintf(fname_interpost_s_2,"interpost_s_%d_%dsec_%d_%d_%d.txt",2,sec,EI_1,EI_2,EI_3);
			sprintf(fname_interpost_s_3,"interpost_s_%d_%dsec_%d_%d_%d.txt",3,sec,EI_1,EI_2,EI_3);
			sprintf(fname_interpost_s_4,"interpost_s_%d_%dsec_%d_%d_%d.txt",4,sec,EI_1,EI_2,EI_3);
			sprintf(fname_interpost_s_5,"interpost_s_%d_%dsec_%d_%d_%d.txt",5,sec,EI_1,EI_2,EI_3);
			sprintf(fname_interpost_s_6,"interpost_s_%d_%dsec_%d_%d_%d.txt",6,sec,EI_1,EI_2,EI_3);


			sprintf(fname_s_1,"s_%d_%dsec_%d_%d_%d.txt",1,sec,EI_1,EI_2,EI_3);
			sprintf(fname_s_2,"s_%d_%dsec_%d_%d_%d.txt",2,sec,EI_1,EI_2,EI_3);
 			sprintf(fname_s_3,"s_%d_%dsec_%d_%d_%d.txt",3,sec,EI_1,EI_2,EI_3);


		  	fp_interpost_s_1=fopen(fname_interpost_s_1,"w");
		  	fp_interpost_s_2=fopen(fname_interpost_s_2,"w");
		  	fp_interpost_s_3=fopen(fname_interpost_s_3,"w");
		  	fp_interpost_s_4=fopen(fname_interpost_s_4,"w");
		  	fp_interpost_s_5=fopen(fname_interpost_s_5,"w");
		  	fp_interpost_s_6=fopen(fname_interpost_s_6,"w");
		  	
			fp_s_1=fopen(fname_s_1,"w");
		  	fp_s_2=fopen(fname_s_2,"w");
			fp_s_3=fopen(fname_s_3,"w");

			
			 for(i=0;i<Ne[0][num];i++)
			  for(j=0;j<interM;j++){
				fprintf(fp_interpost_s_1,"%f\t%d\t%d\n",interpost_s[0][0][i][j],i,interpost[0][0][i][j]);
			       }
			for(i=0;i<Ne[1][num];i++)
			  for(j=0;j<interM;j++){
				fprintf(fp_interpost_s_2,"%f\t%d\t%d\n",interpost_s[1][1][i][j],i,interpost[1][1][i][j]);
			 
			  }
			for(i=0;i<Ne[1][num];i++)
			  for(j=0;j<interM;j++){
				fprintf(fp_interpost_s_3,"%f\t%d\t%d\n",interpost_s[1][0][i][j],i,interpost[1][0][i][j]);
			       }
			for(i=0;i<Ne[2][num];i++)
			  for(j=0;j<interM;j++){
				fprintf(fp_interpost_s_4,"%f\t%d\t%d\n",interpost_s[2][1][i][j]
,i,interpost[2][1][i][j]);
			 
			  }
			for(i=0;i<Ne[2][num];i++)
			  for(j=0;j<interM;j++){
				fprintf(fp_interpost_s_5,"%f\t%d\t%d\n",interpost_s[2][0][i][j]
,i,interpost[2][0][i][j]);
			       }
			for(i=0;i<Ne[0][num];i++)
			  for(j=0;j<interM;j++){
				fprintf(fp_interpost_s_6,"%f\t%d\t%d\n",interpost_s[0][1][i][j]
,i,interpost[0][1][i][j]);
			 
			  }



	                for(i=0;i<Ne[0][num];i++)
			  for(j=0;j<M;j++){
				fprintf(fp_s_1,"%f\t%d\t%d\n",s[/*gr=*/0][i][j],i,post[0][i][j]);	//クッソ長いファイル注意!
			  }
			
		       for(i=0;i<Ne[1][num];i++)
			  for(j=0;j<M;j++){
				fprintf(fp_s_2,"%f\t%d\t%d\n",s[/*gr=*/1/*=Ngr-1*/][i][j],i,post[0][i][j]);	//クッソ長いファイル注意!
			  }

			for(i=0;i<Ne[2][num];i++)
			  for(j=0;j<M;j++){
				fprintf(fp_s_3,"%f\t%d\t%d\n",s[/*gr=*/2/*=Ngr-1*/][i][j],i,post[0][i][j]);	//クッソ長いファイル注意!
			  }

			fclose(fp_interpost_s_1);
 
			fclose(fp_interpost_s_2);

			fclose(fp_interpost_s_3);
 
			fclose(fp_interpost_s_4);

			fclose(fp_interpost_s_5);
 
			fclose(fp_interpost_s_6);

			
			fclose(fp_s_1);
			fclose(fp_s_2);
			fclose(fp_s_3);


		}//for if
//////////////////////////グループ間結合重みヒストグラム作成/////////////////////////////////////////////////     グループ２の重みを調べる
 	//	if(sec==0||sec==100||sec==1000||sec==1999)
	//	{ 
	//		sprintf(fname_interpost_s_1,"interpost_s_%d_%dsec_%d_%d_%d.txt",1,sec,EI_1,EI_2,EI_3);
	//		sprintf(fname_interpost_s_2,"interpost_s_%d_%dsec_%d_%d_%d.txt",2,sec,EI_1,EI_2,EI_3);
	//		sprintf(fname_interpost_s_3,"interpost_s_%d_%dsec_%d_%d_%d.txt",3,sec,EI_1,EI_2,EI_3);
	//		sprintf(fname_interpost_s_4,"interpost_s_%d_%dsec_%d_%d_%d.txt",4,sec,EI_1,EI_2,EI_3);
	//		sprintf(fname_interpost_s_5,"interpost_s_%d_%dsec_%d_%d_%d.txt",5,sec,EI_1,EI_2,EI_3);
	//		sprintf(fname_interpost_s_6,"interpost_s_%d_%dsec_%d_%d_%d.txt",6,sec,EI_1,EI_2,EI_3);


	//		sprintf(fname_s_1,"s_%d_%dsec_%d_%d_%d.txt",1,sec,EI_1,EI_2,EI_3);
	//		sprintf(fname_s_2,"s_%d_%dsec_%d_%d_%d.txt",2,sec,EI_1,EI_2,EI_3);
 	//		sprintf(fname_s_3,"s_%d_%dsec_%d_%d_%d.txt",3,sec,EI_1,EI_2,EI_3);

//
//		  fp_interpost_s_1=fopen(fname_interpost_s_1,"w");
//		  fp_interpost_s_2=fopen(fname_interpost_s_2,"w");
//		  fp_interpost_s_3=fopen(fname_interpost_s_3,"w");
//		  fp_interpost_s_4=fopen(fname_interpost_s_4,"w");
//		  fp_interpost_s_5=fopen(fname_interpost_s_5,"w");
//		  fp_interpost_s_6=fopen(fname_interpost_s_6,"w");



/*		  fp_s_1=fopen(fname_s_1,"w");
		  fp_s_2=fopen(fname_s_2,"w");
		  fp_s_3=fopen(fname_s_3,"w");


		  for(k=0;k<101;k++)interpost_s_1_range[k]=0;		//重みカウント変数の初期化
 		  for(k=0;k<101;k++)interpost_s_2_range[k]=0;		//重みカウント変数の初期化
		  for(k=0;k<101;k++)interpost_s_3_range[k]=0;		//重みカウント変数の初期化
 		  for(k=0;k<101;k++)interpost_s_4_range[k]=0;		//重みカウント変数の初期化
	          for(k=0;k<101;k++)interpost_s_5_range[k]=0;		//重みカウント変数の初期化
 		  for(k=0;k<101;k++)interpost_s_6_range[k]=0;		//重みカウント変数の初期化

 		 for(gr=0;gr<Ngr;gr++)
		  for(k=0;k<101;k++)s_range[gr][k]=0;

		  for(i=0;i<Ne[0][num];i++)
		  {
		    //for(agr=0;agr<Ngr-1;agr++)
                    for(j=0;j<interM;j++)
                    {
		      for(l=0;l<101;l++)
                      {
                         if(l<=interpost_s[0][/*agr==*//*0][i][j]*10 && interpost_s[0][0][i][j]*10<l+1)interpost_s_1_range[l]++;*/
	//	       }
          //          }
            //      }  

		//  for(i=0;i<Ne[1][num];i++)
		  //{
		    //for(agr=0;agr<Ngr-1;agr++)
                    //for(j=0;j<interM;j++)
                    //{
		      //for(l=0;l<101;l++)
                      //{
                        // if(l<=interpost_s[1][1][i][j]*10&&interpost_s[1][1][i][j]*10<l+1)interpost_s_2_range[l]++;
		       //}
                    //}
                 /* }   

		for(i=0;i<Ne[1][num];i++)
		  {
		    //for(agr=0;agr<Ngr-1;agr++)
                    for(j=0;j<interM;j++)
                    {
		      for(l=0;l<101;l++)
                      {
                         if(l<=interpost_s[1][0][i][j]*10&&interpost_s[1][0][i][j]*10<l+1)interpost_s_3_range[l]++;
		       }
                    }
                  }  

		for(i=0;i<Ne[2][num];i++)
		  {
                    for(j=0;j<interM;j++)
                    {
		      for(l=0;l<101;l++)
                      {
                         if(l<=interpost_s[2][1][i][j]*10&&interpost_s[2][1][i][j]*10<l+1)interpost_s_4_range[l]++;
		       }
                    }
                  }  

		for(i=0;i<Ne[0][num];i++)
		  {
                    for(j=0;j<interM;j++)
                    {
		      for(l=0;l<101;l++)
                      {
                         if(l<=interpost_s[0][1][i][j]*10&&interpost_s[0][1][i][j]*10<l+1)interpost_s_5_range[l]++;
		       }
                    }
                  }  

		for(i=0;i<Ne[2][num];i++)
		  {
                    for(j=0;j<interM;j++)
                    {
		      for(l=0;l<101;l++)
                      {
                         if(l<=interpost_s[2][0][i][j]*10&&interpost_s[2][0][i][j]*10<l+1)interpost_s_6_range[l]++;
		       }
                    }
                  }  


//////////////////////////////////                 ///////////////////////////////////////////////////////////////////////////////

                    for(gr=0;gr<Ngr;gr++)
                    for(i=0;i<Ne[gr][num];i++)
		  {
                    for(j=0;j<M;j++)
                    {
		      for(l=0;l<101;l++)
                      {
                         if(l<=s[gr][i][j]*10&&s[gr][i][j]*10<l+1)s_range[gr][l]++;
		       }
                    }
                  }      

//////////////////////////////////    正規化したヒストグラムを出力する	///////////////////////////////////////////////////////////
			for(l=0;l<101;l++)
                        {
                         interpost_s_1_pd[l]=float(interpost_s_1_range[l])/float(Ne[0][num])/float(interM);
                         fprintf(fp_interpost_s_1,"%d\t%f\n",l,interpost_s_1_pd[l]);		//割合(確率密度)を出力

			}
                   fclose(fp_interpost_s_1);
			for(l=0;l<101;l++)
                        {
                         interpost_s_2_pd[l]=float(interpost_s_2_range[l])/float(Ne[1][num])/float(interM);
                         fprintf(fp_interpost_s_2,"%d\t%f\n",l,interpost_s_2_pd[l]);		//割合(確率密度)を出力

			}


                   fclose(fp_interpost_s_2);

			for(l=0;l<101;l++)
                        {
                         interpost_s_3_pd[l]=float(interpost_s_3_range[l])/float(Ne[1][num])/float(interM);
                         fprintf(fp_interpost_s_3,"%d\t%f\n",l,interpost_s_3_pd[l]);		//割合(確率密度)を出力

			}
                   fclose(fp_interpost_s_3);
			for(l=0;l<101;l++)
                        {
                         interpost_s_4_pd[l]=float(interpost_s_4_range[l])/float(Ne[2][num])/float(interM);
                         fprintf(fp_interpost_s_4,"%d\t%f\n",l,interpost_s_4_pd[l]);		//割合(確率密度)を出力

			}


                   fclose(fp_interpost_s_4);
			for(l=0;l<101;l++)
                        {
                         interpost_s_5_pd[l]=float(interpost_s_5_range[l])/float(Ne[0][num])/float(interM);
                         fprintf(fp_interpost_s_5,"%d\t%f\n",l,interpost_s_5_pd[l]);		//割合(確率密度)を出力

			}
                   fclose(fp_interpost_s_5);

			for(l=0;l<101;l++)
                        {
                         interpost_s_6_pd[l]=float(interpost_s_6_range[l])/float(Ne[2][num])/float(interM);
                         fprintf(fp_interpost_s_6,"%d\t%f\n",l,interpost_s_6_pd[l]);		//割合(確率密度)を出力

			}
                   fclose(fp_interpost_s_6);




			
			for(i=0;i<101;i++)
                        {
                         s_pd[0][i]=float(s_range[0][i])/float(Ne[0][num])/float(M);
                         fprintf(fp_s_1,"%d\t%f\n",i,s_pd[0][i]);		//割合(確率密度)を出力
			}

                   fclose(fp_s_1);

 			for(i=0;i<101;i++)
                        {
                         s_pd[1][i]=float(s_range[1][i])/float(Ne[1][num])/float(M);
                         fprintf(fp_s_2,"%d\t%f\n",i,s_pd[1][i]);		//割合(確率密度)を出力
			}

                   fclose(fp_s_2);

			for(i=0;i<101;i++)
                        {
                         s_pd[2][i]=float(s_range[2][i])/float(Ne[2][num])/float(M);
                         fprintf(fp_s_3,"%d\t%f\n",i,s_pd[2][i]);		//割合(確率密度)を出力
			}

                   fclose(fp_s_3);


                 }//ifのfor文終わり
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////

   	}//secについてのfor文終わり　//for14　end




 
	              
        

       for(i=0;i<Ngr;i++){
       printf("%d %d\n",Ne[i][num],Ni[i][num]);
       
        }


	



	//MI_1=mutual_information(LAPx, LAPy,1000);
//	MI_2=mutual_information(LAPy, LAPz,1000);
//	MI_3=mutual_information(LAPz, LAPx,1000);


//	fp_MI=fopen("MI_index.txt","a");
//	fprintf(fp_MI,"%d\t%d\t%d\t\t%f\t%f\t%f\n",EI_1,EI_2,EI_3,MI_1,MI_2,MI_3);
//	fclose(fp_MI);





	 if(num<NUM_EI_PATTERN-1)cout <<"次のEIの組み合わせのシミュレーションを実行します"<<endl;
	 else	    cout <<"シミュレーション終了です"<<endl;

    }//numの終わり
      
 fp_firing_count=fopen("firing_count.txt","w");
        for(gr=0;gr<Ngr;gr++)
         for(i=0;i<N;i++)
          for(num=0;num<NUM_EI_PATTERN;num++)
         
          if(num<NUM_EI_PATTERN-1)fprintf(fp_firing_count,"%d\t",firing_count[num][gr][i]);
          else fprintf(fp_firing_count,"%d\n",firing_count[num][gr][i]);
          if(i==N-1)fprintf(fp_firing_count,"\n\n");

        fclose(fp_firing_count);

       fp_firing_rate=fopen("firing_rate_seq.txt","w");
         for(gr=0;gr<Ngr;gr++)
          for(sec=0;sec<SIM_TIME;sec++)
           for(num=0;num<NUM_EI_PATTERN;num++)
            if(num<NUM_EI_PATTERN-1)fprintf(fp_firing_rate,"%f\t",firing_rate_seq[num][gr][sec]);
            else fprintf(fp_firing_rate,"%f\n",firing_rate_seq[num][gr][sec]);
            

               fclose(fp_firing_rate);

        system("mv spikes* data_spikes"); 
	system("mv LAP* data_LAP");
	system("mv s_* data_s");
	system("mv interpost_s* data_interpost_s");
	 
        system("mv firing_count.txt data_firing_count");
        system("mv firing_rate_seq.txt data_firing_count");
	system(str); 

	
 
     
     }//trial文  
       t2=time(NULL);	//終了時間
       printf("time=%d[s]\n",(int)(t2-t1));
}	//main関数の終わり
