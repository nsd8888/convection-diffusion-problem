#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define grid 50

using namespace std;
class simple
{

private:
public:
    double a_E, a_P, a_N, a_S;
    double ug[grid][grid + 1],un[grid][grid + 1], ustr[grid][grid + 1], u_E, u_P, u_N, u_S, a_e[grid][grid+1], uc[grid][grid];

    double vg[grid + 1][grid],vn[grid + 1][grid], vstr[grid + 1][grid], v_E, v_W, v_N, v_S, b_e[grid+1][grid], vc[grid][grid];
    double b_E, b_W, b_N, b_S;

    double p[grid + 1][grid + 1], pn[grid + 1][grid + 1], pstr[grid + 1][grid + 1], p_e[grid + 1][grid + 1], pcor[grid + 1][grid + 1], c_p[grid + 1][grid + 1], pc[grid][grid];

    int i, j, a;
    double Re = 100.0;
    double alpha = 0.8;
    double dx = 1.0 / (grid - 1);
    double dy = 1.0 / (grid - 1);

    //error = 0.0;
    void initialise();
    void operate();
    void Pr_correction();
    void Pr_BC();
    void vel_BC();
    void calcerror();
    void iterate();
};

void simple::initialise()
{
    for (i = 0; i <= (grid); i++)
    {
        for (j = 0; j <= (grid); j++)
        {
           //pstr[i][j] = 1 + ((9810 * ((grid - 1) - ((j - 1) + (1 / 2)) * dy) * 0.00001));
           pstr[i][j]=0.5;
            //pstr[0][j] = pstr[1][j];
           // pstr[grid][j] = pstr[grid - 1][j];
        }
        //pstr[i][grid] = pstr[i][grid - 1];
    }

    //initialize the velocity field
    for (i = 0; i <= (grid - 1); i++)
    {
        
        for (j = 0; j <= (grid - 2); j++)
        {
            ug[i][j] = 0.0;
        }
        ug[i][grid] = 1.0;
        ug[i][grid - 1] = 1.0;
    }

    //cout<<ug[0][2];

    for (i = 0; i <= (grid); i++)
    {
        for (j = 0; j <= (grid - 1); j++)
        {
            vg[i][j] = 0.0;
        }
    }
  //cout<<vg[2][1];
    cout <<"variable initialization is completed"<<endl;
}

void simple::operate()
{

    for (i =(grid-2); i >= 1; i--)
    {
        for (j = (grid-1); j >= 1; j--)
        {

            u_E = (ug[i + 1][j] + ug[i][j]) / 2;
            u_P = (ustr[i][j] + ug[i - 1][j]) / 2;
            u_N = (vg[i][j] + vg[i + 1][j]) / 2;
            u_S = (vg[i][j - 1] + vg[i + 1][j - 1]) / 2;

            a_E = ((dy / dx * Re) - (u_E * dy / 2));
            a_P = ((dy / Re * dx) + (u_P * dy / 2));
            a_N = ((dx / dy * Re) - (u_N * dx / 2));
            a_S = ((dx / Re * dy) + (u_S * dx / 2));

            a_e[i][j] = ((0.5 * u_E) - (0.5 * u_P)) * dy + ((0.5 * u_N) - (0.5 * u_S)) * dx + ((dy / dx) * (2 / Re)) + ((dx / dy) * (2 / Re));

            ustr[i][j] = ((a_E / (a_e[i][j]+1.0) * ug[i + 1][j])) + ((a_P / (a_e[i][j]+1.0) * ug[i - 1][j])) + ((a_N / (a_e[i][j]+1.0) * ug[i][j + 1]))+ ((a_S / (a_e[i][j]+1.0) * ug[i][j - 1]))
                         + ((dy / (a_e[i][j]+1.0) * (pstr[i + 1][j] - pstr[i][j])));

        }
        
     }

        for (j = 1; j <= (grid - 1); j++)
        {
            ustr[0][j] = 0;
            ustr[grid][j] = 0;
        }

        for (i = 0; i <= (grid - 1); i++)
        {
            ustr[i][0] = -ustr[i][1];
            ustr[i][grid] = 2 - ustr[i][grid - 1];
        }


     for(i = (grid-1); i >= 1; i--)
     {
        for (j =(grid-2); j >= 1; j--)
        {
            v_E = (ug[i][j + 1] + ug[i][j]) / 2;
            v_W = (ug[i - 1][j] + ug[i - 1][j + 1]) / 2;
            v_N = (vg[i][j + 1] + vg[i][j]) / 2;
            v_S = (vg[i][j] + vg[i][j - 1]) / 2;

            b_W = ((v_W * (dy / 2)) + (dy / (dx * Re)));
            b_E = ((dy / (dx * Re)) - (v_E * (dy / 2)));
            b_N = ((dx / (Re * dy)) - (v_N * (dx / 2)));
            b_S = ((dx / (Re * dy)) + (v_S * (dx / 2)));

            b_e[i][j] = ((v_E - v_W) * (dy / 2)) + ((v_N - v_S) * (dx / 2)) + ((2 / Re) * (dy / dx)) + ((2 / Re) * (dx / dy));

            vstr[i][j] = (b_W / (b_e[i][j]+1.0) * vg[i - 1][j]) + (b_E / (b_e[i][j]+1.0) * vg[i + 1][j]) + (b_N / (b_e[i][j]+1.0) * vg[i][j + 1] )+ (b_S / (b_e[i][j]+1.0) * vg[i][j - 1]) + (pstr[i][j + 1] - pstr[i][j]) * (dx / (b_e[i][j]+1.0));
        
          
        }
    }
     // cout<<b_e[1][1]<<endl;
      //cout<<vstr[1][1]<<endl;
    for (j = 1; j <= (grid - 2); j++)
    {
        vstr[0][j] = -vstr[1][j];
        vstr[grid][j] = -vstr[grid - 1][j];
    }

    for (i = 0; i <= (grid); i++)
    {
        vstr[i][0] = 0.0;
        vstr[i][grid - 1] = 0.0;
    }
 // cout<<a_e[2][1]<<endl;
  //cout<<a_e[3][1]<<endl;
  //cout<<ustr[2][1]<<endl;
  //cout<<ustr[3][1]<<endl;

}

void simple::Pr_correction()
{
    //pressure correction

    //first initialie them as the zero
    for (i = 0; i <= (grid); i++)
    {
        for (j = 0; j <= (grid); j++)
        {
            pcor[i][j] = 0;
        }
    }

    //calculate the pressure corection at interior points
    for (i =(grid-1); i >= 1; i--)
    {
        for (j = (grid-1); j >= 1; j--)
        {
            p_e[i][j] = -(dy * dy / (a_e[i - 1][j]+1.0)) - (dy * dy / (a_e[i][j]+1.0)) - (dx * dx / (b_e[i][j]+1.0)) - (dx * dx / (b_e[i][j - 1]+1.0));
            c_p[i][j] = (dy * (ustr[i - 1][j] - ustr[i][j])) + (dx * (vstr[i][j - 1] - vstr[i][j]));
            pcor[i][j] = (-dy * dy / ((a_e[i][j]+1.0) * p_e[i][j])) * pcor[i + 1][j] + (-dy * dy / ((a_e[i - 1][j]+1.0) * p_e[i][j]) * pcor[i - 1][j]) + (-dx * dx / ((b_e[i][j]+1.0) * p_e[i][j]) * pcor[i][j + 1]) + (-dx * dx / ((b_e[i][j - 1]+1.0) * p_e[i][j]) * pcor[i][j - 1]) + c_p[i][j];
        }
    }
    //cout<<p_e[48][48]<<endl;
    //cout<<c_p[48][48]<<endl;
    //cout<<pcor[48][48]<<endl;
}

void simple::Pr_BC()
{

    //interior pressure boundaru condition

    for (i = 1; i <= (grid - 1); i++)
    {
        for (j = 1; j <= (grid - 1); j++)
        {
            pn[i][j] = pstr[i][j] + (alpha * pcor[i][j]);
        }
    }

    //boundary p
    for (i = 1; i <= (grid - 1); i++)
    {
        pn[i][0] = pn[i][1];
        pn[i][grid] = pn[i][grid - 1];
    }

    for (j = 0; j <= (grid); j++)
    {
        pn[0][j] = pn[1][j];
        pn[grid][j] = pn[grid - 1][j];
    }


}

void simple::vel_BC()
{
    //update u and v values
    for (i = 1; i <= (grid - 2); i++)
    {
        for (j = 1; j <= (grid - 1); j++)
        {
            un[i][j] = ustr[i][j] + (alpha * (dy/(a_e[i][j]+1.0)) * (pcor[i + 1][j] - pcor[i][j]));
        }
    }

    for (i = 1; i <= (grid - 1); i++)
    {
        for (j = 1; j <= (grid - 2); j++)
        {
            vn[i][j] = vstr[i][j] + (alpha * ((dy/b_e[i][j]+1.0)) * (pcor[i][j + 1] - pcor[i][j]));
        }
    }
    //update the boundary conditions

    for (j = 1; j <= (grid - 1); j++)
    {
        un[0][j] = 0;
        un[grid][j] = 0;
    }

    for (i = 0; i <= (grid - 1); i++)
    {
        un[i][0] = -un[i][1];
        un[i][grid] = 2 - un[i][grid - 1];
    }

    for (j = 1; j <= (grid - 2); j++)
    {
        vn[0][j] = -vn[1][j];
        vn[grid][j] = -vn[grid - 1][j];
    }

    for (i = 0; i <= (grid); i++)
    {
        vn[i][0] = 0.0;
        vn[i][grid - 1] = 0.0;
    }
}

void simple::calcerror()
{

    //error
    double error = 0.0;
    for (i = 1; i <= (grid - 1); i++)
    {
        for (j = 1; j <= (grid - 1); j++)
        {
            error = error + abs(c_p[i][j]);
        }
    }

    cout << "the error at the step is " << error << endl;
}

void simple::iterate()
{
    // Iterating u
    for (i = 0; i <= (grid - 1); i++)
    {
        for (j = 0; j <= (grid); j++)
        {
            ug[i][j] = un[i][j];
        }
    }

    // Iterating v
    for (i = 0; i <= (grid); i++)
    {
        for (j = 0; j <= (grid - 1); j++)
        {
            vg[i][j] = vn[i][j];
        }
    }

    // Iterating p
    for (i = 0; i <= (grid); i++)
    {
        for (j = 0; j <= (grid); j++)
        {
            pstr[i][j] = pn[i][j];
        }
    }
}

int main()
{
    simple lid;
    int s, n;
    lid.initialise();
    
    cout<<"enter the value for number of steps"<<endl;
    cin>>s;
    for(n=1; n<=s; n++)
    {
    lid.operate();
    lid.Pr_correction();
    lid.Pr_BC();
    lid.vel_BC();
    lid.calcerror();
    lid.iterate();
    }
    return 0;
}
