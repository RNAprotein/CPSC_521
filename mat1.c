#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "mpi.h"

//main函数代码：
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR



int master();
int slave();

int main(void) {
	
    int test;
	int myid; //
	int numprocs; //	
	MPI_Status status;
	
	//MPI程序的初始化
	MPI_Init(NULL,NULL);
	//获取当前进程的进程编号，保存在myid中
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	//获取当前进程的进程总数，保存在numprocs中
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	
	//如果进程数小于2个，则退出程序
	if(numprocs<2)
	{			
		cout<<"Too Few Processes!"<<endl;		
		MPI_Abort(MPI_COMM_WORLD,99);
	}
	
	if(myid==0)
	{
		//id是0的进程作为主进程
		master();
	}
	else
	{
		//id是非0的进程作为从进程
		slave();
	}
	
	//MPI程序释放资源
	MPI_Finalize();
	return EXIT_SUCCESS;
}


//主进程
int master()
{
	int myid; //
	int numprocs; //	
	MPI_Status status;
	
	//获取当前进程的进程编号，保存在myid中
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	//获取当前进程的进程总数，保存在numprocs中
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

	const int M=4;
        const int K=3;
        const int N=2;
	//结果矩阵C
	int C[M*N];
	int V[M];
	//将结果矩阵各元素初始化为0
	for(int i=0;i<M*N;i++)
	{
		C[i]=0;
	}
	//C矩阵共N列数据
	for(int i=0;i<N;i++)
	{		
		//从各从进程接受数据
		for(int j=1;j<numprocs;j++)
		{
			MPI_Recv(V, M, MPI_INT, j, i, MPI_COMM_WORLD, &status);	
			//将各从进程的部分数据整合成完整的一个向量
			for(int ij=0;ij<M;ij++)
			{
				C[ij*N+i]=C[ij*N+i]+V[ij];
			}
		}		
	}
	//打印结果
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)
		{
			cout<<C[i*N+j]<<"\t";
		}
		cout<<endl;
	}	
	return 0;
}
//从进程代码：
//从进程
int slave()
{
	int myid; //
	int numprocs; //	
	MPI_Status status;
	
	//获取当前进程的进程编号，保存在myid中
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	//获取当前进程的进程总数，保存在numprocs中
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	
	int const M(4),K(3),N(2);
	vector<int> indexes;	

	int V[M];
	
	int A[M*K]={1,2,3,4,5,6,7,8,9,8,7,6};
	int B[K*N]={5,4,3,2,1,0};
	 
	//计算本进程分配A矩阵中哪些行，保存在indexes中
	//比如2个从进程，A是4行数据，id为1的从进程分配的是0和2行
	for(int i=myid-1;i<M;i+=numprocs-1)
	{
		indexes.push_back(i);
	} 
	//依次计算B的N列数据
	for(int i=0;i<N;i++)
	{
		//将保存列数据的向量V各元素初始化为0
		for(int j=0;j<M;j++)
		{
			V[j]=0;
		}
		if(indexes.size()>0)
		{						
			for(int j=0;j<indexes.size();j++)
			{
				for(int ij=0;ij<K;ij++)
				{
					V[indexes[j]]+=A[indexes[j]*K+ij]*B[ij*N+i];
				}
			}
			MPI_Send(V, M, MPI_INT, 0,i,MPI_COMM_WORLD);				
		}
		else
		{
			//发送全0的V
			MPI_Send(V, M, MPI_INT, 0,i,MPI_COMM_WORLD);	
		}
	}	
	return 0;
}

