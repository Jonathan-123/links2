#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define pi 3.14159265359

/*
	To Do:
		fb -done
		cs -done
		ics -done
		pAS -NA
		pBP	-NA
		pBU	-NA
		tran -NA
		grash -NA
		extran -NA
		tog -NA
		fbv -done
		csv -done
		csa	-done
		fba -done
*/

double PI(){
	return pi;
}

void  fb(float a, float b, float c, float d, float T1, float T2, float *T3, float *T4, bool crossed){
		float t2 = (T2-T1)*pi/180;
		float K1 = d/a;
		float K2 = d/c;
		float K3 = (pow(a,2)-pow(b,2)+pow(c,2)+pow(d,2))/(2*a*c);
		float K4 = d/b;
		float K5 = (pow(c,2)-pow(d,2)-pow(a,2)-pow(b,2))/(2*a*b);
		float A = cos(t2)-K1-K2*cos(t2)+K3;
	    float B = -2*sin(t2);
	    float C = K1-(K2+1)*cos(t2)+K3;
		float D = cos(t2)-K1+K4*cos(t2)+K5;
		float E = -2*sin(t2);
		float F = K1+(K4-1)*cos(t2)+K5;
		
		float t3=0, t4=0;
		//Normal Configuration
		if (crossed == false){
			t3 = 2*atan((-E-sqrt(pow(E,2)-4*D*F))/(2*D));
			t4 = 2*atan((-B-sqrt(pow(B,2)-4*A*C))/(2*A));
		}
		//Crossed Configuration
		if (crossed == true){
			t3 = 2*atan((-E+sqrt(pow(E,2)-4*D*F))/(2*D));
        	t4 = 2*atan((-B+sqrt(pow(B,2)-4*A*C))/(2*A));

		}
		
		*T3 = t3*180/pi+T1;
		*T4 = t4*180/pi+T1;
		//*T31 = t31*180/pi+T1;
		//*T32 = t32*180/pi+T1;
		//*T41 = t41*180/pi+T1;
		//*T42 = t42*180/pi+T1;
}
void cs(float a, float b, float c, float T1, float T2, float *T3, float *d, bool crossed){
		float t2 = (T2-T1)*pi/180;
		if (crossed == false){
			float t3 = asin((a*sin(t2)-c)/b);
			*d = a*cos(t2) - b*cos(t3);
			*T3 = t3*180/pi+T1;
		}
		if (crossed == true){
			float t3 = pi - asin((a*sin(t2)-c)/b);
			*d = a*cos(t2) - b*cos(t3);
			*T3 = t3*180/pi+T1;
		}
}
void ics(float a, float c, float d, float T2, float G, float *T3, float *T4, float *b, bool crossed){
		float t2 = T2*pi/180;
		float g = G*pi/180;
		float P = a*sin(t2)*sin(g)+(a*cos(t2)-d)*cos(g);
		float Q = -a*sin(t2)*cos(g)+(a*cos(t2)-d)*sin(g);
		float R = -c*sin(g);
		float S = R-Q;
		float T = 2*P;
		float U = Q+R;
		float t3,t4;
		if (crossed == false){
			t4 = 2*atan((-T-sqrt(pow(T,2)-4*S*U))/(2*S));
			t3 = pi+t4+g;
			*b = (a*sin(t2)-c*sin(t4))/sin(t3);
			*T3 = 180*t3/pi;
			*T4 = 180*t4/pi;
		}
		if (crossed == true){
			t4 = 2*atan((-T+sqrt(pow(T,2)-4*S*U))/(2*S));
			t3 = t4+g;
			*b = (a*sin(t2)-c*sin(t4))/sin(t3);
			*T3 = 180*t3/pi;
			*T4 = 180*t4/pi;
		}	
}

/*void fbpos(float a, float b, float c, float d, float s, float p, float u, float T1, float T2, float T3, float T4, float O2x, float O2y, float &Ax, float &Ay, float &Bx, float &By, float &Sx, float
*/
void pAS(float O2x, float O2y, float T2, float d2, float s, float *Sx, float *Sy){
		*Sx = O2x+s*cos((T2+d2)*pi/180);
		*Sy = O2y+s*sin((T2+d2)*pi/180);
}
void pBP(float Ax, float Ay, float T3, float d3, float p, float *Px, float *Py){
		*Px = Ax+p*cos((T3+d3)*pi/180);
		*Py = Ay+p*sin((T3+d3)*pi/180);
}
void pBU(float O4x, float O4y, float T4, float d4, float u, float *Ux, float *Uy){
		*Ux = O4x+u*cos((T4+d4)*pi/180);
		*Uy = O4y+u*sin((T4+d4)*pi/180);
}
void tran(float T3, float T4, float *mu){
		if(sqrt(pow(T4-T3,2))<90)
		{
			 *mu = sqrt(pow(T4-T3,2));
			 }
		else{
			 *mu = 180 - sqrt(pow(T4-T3,2));
			 }
}
int grash(float a, float b, float c, float d){
		//float lengths[4] = {a,b,c,d};
		float l = fmaxf(fmaxf(a,b),fmaxf(c,d));
		float s = fminf(fminf(a,b),fminf(c,d));
		float p = fmaxf(fminf(a,b),fminf(c,d));
		float q = fminf(fmaxf(a,b),fmaxf(c,d));
		int x = 0;
		if(s+l<=p+q){
			x = 1;
		}
		return x;
}
void extran(float a, float b, float c, float d, float *mu1, float *mu2){
	if(grash(a,b,c,d)){	
		*mu1 = acos((pow(b,2)+pow(c,2)-pow(d-a,2))/(2*b*c))*180/pi;
		*mu2 = 180 - acos((pow(b,2)+pow(c,2)-pow(d+a,2))/(2*b*c))*180/pi;
}	
	else{
		*mu1 = 0;
		*mu2 = acos((pow(a+b,2)+pow(c,2)-pow(d,2))/(2*c*(a+b)))*180/pi;
		}
}
void tog(float a, float b, float c, float d, float *Ttog1, float *Ttog2){
	if(!grash(a,b,c,d)){		
		float s = (pow(a,2)+pow(d,2)-pow(b,2)-pow(c,2))/(2*a*d)+(b*c)/(a*d);
		float t = (pow(a,2)+pow(d,2)-pow(b,2)-pow(c,2))/(2*a*d)-(b*c)/(a*d);
		if(pow(s,2)<=1){
				*Ttog1 = acos(s)*180/pi;
				*Ttog2 = -1*acos(s)*180/pi;
		}
		else{
				*Ttog1 = acos(t)*180/pi;
				*Ttog2 = -1*acos(t)*180/pi;
		}
	}
	else{
			printf("The linkage is grashof\n");
	}
}
void fbv(float a, float b, float c, float T2, float T3, float T4, float w2, float *w3, float *w4, bool crossed){
	T2 = T2*pi/180;
	if (crossed == false){
		T3 = T3*pi/180;
		T4 = T4*pi/180;
		*w3 = (a*w2/b)*sin(T4-T2)/sin(T3-T4);
		*w4 = (a*w2/c)*sin(T2-T3)/sin(T4-T3);
	}
	if (crossed == true){
		T3 = T3*pi/180;
		T4 = T4*pi/180;
		*w3 = (a*w2/b)*sin(T4-T2)/sin(T3-T4);
		*w4 = (a*w2/c)*sin(T2-T3)/sin(T4-T3);
	}
}
void csv(float a, float b, float c, float T1, float T2, float T3, float w2, float *w3, float *ddt, bool crossed){
	T1 = T1*pi/180;
	T2 = T2*pi/180;
    if(crossed ==false){
		T3 = T3*pi/180;
		*w3 = (a/b)*(cos(T2-T1)/cos(T3-T1))*w2;
		*ddt = -a*w2*sin(T2-T1)+b**w3*sin(T3-T1);

	}
	if(crossed == true){
		T3 = T3*pi/180;
		*w3 = (a/b)*(cos(T2-T1)/cos(T3-T1))*w2;
		*ddt = -a*w2*sin(T2-T1)+b**w3*sin(T3-T1);
	}
	
    
}
void csa(float a, float b, float c, float T1, float T2, float T3, float w2, float w3, float al2, float *al3, float *ddtt, bool crossed){
	T1 = T1*pi/180;
	T2 = T2*pi/180;
	if(crossed == false){
		T3 = T3*pi/180;
		*al3 = (a*al2*cos(T2-T1)-a*pow(w2,2)*sin(T2-T1)+b*pow(w3,2)*sin(T3-T1))/(b*cos(T3-T1));
		*ddtt = -a*al2*sin(T2-T1)-a*pow(w2,2)*cos(T2-T1)+b**al3*sin(T3-T1)+b*pow(w3,2)*cos(T3-T1);
	}
	if(crossed == true){
		T3 = T3*pi/180;
		*al3 = (a*al2*cos(T2-T1)-a*pow(w2,2)*sin(T2-T1)+b*pow(w3,2)*sin(T3-T1))/(b*cos(T3-T1));
		*ddtt = -a*al2*sin(T2-T1)-a*pow(w2,2)*cos(T2-T1)+b**al3*sin(T3-T1)+b*pow(w3,2)*cos(T3-T1);
	}
	   
    
    
    
    
    
}


void fba(float a, float b, float c, float T2, float T3, float T4, float w2, float w3, float w4, float al2, float *al3, float *al4, bool crossed){
	T2 = T2*pi/180;
	if(crossed == false){
		T3 = T3*pi/180;
		T4 = T4*pi/180;
		float A1 = c*sin(T4);
		float B1 = b*sin(T3);
		float C1 = a*al2*sin(T2)+a*pow(w2,2)*cos(T2)+b*pow(w3,2)*cos(T3)-c*pow(w4,2)*cos(T4);
		float D1 = c*cos(T4);
		float E1 = b*cos(T3);
		float F1 = a*al2*cos(T2)-a*pow(w2,2)*sin(T2)-b*pow(w3,2)*sin(T3)+c*pow(w4,2)*sin(T4);
		*al3 = (C1*D1-A1*F1)/(A1*E1-B1*D1);
		*al4 = (C1*E1-B1*F1)/(A1*E1-B1*D1);
	}
	if(crossed == true){
		T3 = T3*pi/180;
		T4 = T4*pi/180;
		float A2 = c*sin(T4);
		float B2 = b*sin(T3);
		float C2 = a*al2*sin(T2)+a*pow(w2,2)*cos(T2)+b*pow(w3,2)*cos(T3)-c*pow(w4,2)*cos(T4);
		float D2 = c*cos(T4);
		float E2 = b*cos(T3);
		float F2 = a*al2*cos(T2)-a*pow(w2,2)*sin(T2)-b*pow(w3,2)*sin(T3)+c*pow(w4,2)*sin(T4);
		*al3 = (C2*D2-A2*F2)/(A2*E2-B2*D2);
		*al4 = (C2*E2-B2*F2)/(A2*E2-B2*D2);
	}
	else{
		fprintf(stderr,"crossed needs to be a boolean (either TRUE and FALSE, or 1 and 0)");
	}
}
