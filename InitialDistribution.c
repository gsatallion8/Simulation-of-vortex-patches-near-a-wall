154.1#include <stdio.h>
#include <math.h>
#include <stdlib.h>

main()
{
    int i,j,k,N,N1,N2,count;
    float *points[2],pi,r1,r2,x1,x2,y1,y2,g,dt,nSteps;
    FILE*fp;
    fp=fopen("VortexMerging_Wall+Doublet.txt","w");

//    printf("Enter the number of vortices required to be distributed per patch\n");
//    scanf("%d",&N);
//    printf("Enter the radial discretization required in the distributions\n");
//    scanf("%d",&N1);
//    printf("Enter the radii of the distributions\n");
//    scanf("%f %f",&r1,&r2);
//    printf("Enter the center of the two distributions\n");
//    scanf("%f %f %f %f",&x1,&y1,&x2,&y2);
//    printf("Enter the strength of the point vortices to be created\n");
//    scanf("%f",&g);
//    printf("Enter the time step and the number of time steps\n");
//    scanf("%f %f",&dt,&nSteps);
    //preferrably N is a multiple of N1*(N1+1)

    N=264;
    N1=11;
    r1=1;
    r2=1;
    x1=0;
    x2=2;
    y1=5;
    y2=7;
    g=10;
    dt=0.00001;
    nSteps=100000;

    points[0]=(float*)malloc(2*N*sizeof(float));
    points[1]=(float*)malloc(2*N*sizeof(float));
    //printf("1\n");
    pi=3.141593;
    count=0;
    for(i=1;i<=N1;i++)
    {
        //printf("2\n");
        N2=2*N*i/(N1*(N1+1));
        printf("--%d--\n",N2);
        for(j=0;j<N2;j++)
        {
            points[0][count]=x1+r1*i*cos(2*pi*j/N2)/N1;
            points[1][count]=y1+r1*i*sin(2*pi*j/N2)/N1;
            count+=1;
            printf("%d\n",count);
        }
    }

    for(i=1;i<=N1;i++)
    {
        //printf("3\n");
        N2=2*N*i/(N1*(N1+1));
        for(j=0;j<N2;j++)
        {
            points[0][count]=x2+r2*i*cos(2*pi*j/N2)/N1;
            points[1][count]=y2+r2*i*sin(2*pi*j/N2)/N1;
            count+=1;
        }
    }

    printf("4\n");
    //fprintf(fp,"Hi!!");
    fprintf(fp,"%d\t",2*N);
    fprintf(fp,"%f\t",dt);
    fprintf(fp,"%f\n",nSteps);
    for(i=0;i<count/2;i++)
    {
        printf("-> %d %d\n",i,count);
        fprintf(fp,"%f\t%f\t%f\n",-g,points[0][i],points[1][i]);
    }
    for(i=count/2;i<count;i++)
    {
        printf("-> %d %d\n",i,count);
        fprintf(fp,"%f\t%f\t%f\n",g,points[0][i],points[1][i]);
    }
    //printf("5\n");
    fclose(fp);
    return 0;
}
