namespace Utils{
	using System;
	[Serializable]
	public class Mat3x3 {
		public double[] data=new double[9];
		public double this[int x, int y] { get { return data[x+3*y]; } set { data[x+3*y]=value; } }
		public Mat3x3() { Set0(); }
		public Mat3x3(Mat3x3 m) { for (int i=0; i<9; i++)data[i]=m.data[i]; }
		public Mat3x3(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8) { data[0]=x0; data[1]=x1; data[2]=x2; data[3]=x3; data[4]=x4; data[5]=x5; data[6]=x6; data[7]=x7; data[8]=x8; }
		public void Mul(double k) { for (int i=0; i<9; i++)data[i]*=k; }
		public void Add(Mat3x3 m) { for (int i=0; i<9; i++)data[i]+=m.data[i]; }
		public void AA(double k, Mat3x3 m) { for (int i=0; i<9; i++)data[i]=k*m.data[i]; }
		public void MA(double k, Mat3x3 m) { for (int i=0; i<9; i++)data[i]+=k*m.data[i]; }
		public void Set(Mat3x3 m) { for (int i=0; i<9; i++)data[i]=m.data[i]; }
		public void Set(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8) { data[0]=x0; data[1]=x1; data[2]=x2; data[3]=x3; data[4]=x4; data[5]=x5; data[6]=x6; data[7]=x7; data[8]=x8; }
		public void Set0() { for (int i=0; i<9; i++)data[i]=0; }
		public void GetAxisVectors(out Vec3 Vx, out Vec3 Vy, out Vec3 Vz) { Vx.x=data[0]; Vy.x=data[1]; Vz.x=data[2]; Vx.y=data[3]; Vy.y=data[4]; Vz.y=data[5]; Vx.z=data[6]; Vy.z=data[7]; Vz.z=data[8]; }
		public void Transpose() {
			double tmp=data[1]; data[1]=data[3]; data[3]=tmp;
			tmp=data[2]; data[2]=data[6]; data[6]=tmp;
			tmp=data[5]; data[5]=data[7]; data[7]=tmp;
		}
		public static void Transpose(Mat3x3 m, Mat3x3 ans) {
			ans.data[0]=m.data[0]; ans.data[1]=m.data[3]; ans.data[2]=m.data[6];
			ans.data[3]=m.data[1]; ans.data[4]=m.data[4]; ans.data[5]=m.data[7];
			ans.data[6]=m.data[2]; ans.data[7]=m.data[5]; ans.data[8]=m.data[8];
		}
		public double Det() { return data[0]*(data[4]*data[8]-data[5]*data[7])+data[1]*(data[5]*data[6]-data[3]*data[8])+data[2]*(data[3]*data[7]-data[4]*data[6]); }
		public static void Prod(Mat3x3 A, Mat3x3 B, Mat3x3 Ans) {
			Ans.data[0]=A.data[0]*B.data[0]+A.data[1]*B.data[3]+A.data[2]*B.data[6];
			Ans.data[1]=A.data[0]*B.data[1]+A.data[1]*B.data[4]+A.data[2]*B.data[7];
			Ans.data[2]=A.data[0]*B.data[2]+A.data[1]*B.data[5]+A.data[2]*B.data[8];
			Ans.data[3]=A.data[3]*B.data[0]+A.data[4]*B.data[3]+A.data[5]*B.data[6];
			Ans.data[4]=A.data[3]*B.data[1]+A.data[4]*B.data[4]+A.data[5]*B.data[7];
			Ans.data[5]=A.data[3]*B.data[2]+A.data[4]*B.data[5]+A.data[5]*B.data[8];
			Ans.data[6]=A.data[6]*B.data[0]+A.data[7]*B.data[3]+A.data[8]*B.data[6];
			Ans.data[7]=A.data[6]*B.data[1]+A.data[7]*B.data[4]+A.data[8]*B.data[7];
			Ans.data[8]=A.data[6]*B.data[2]+A.data[7]*B.data[5]+A.data[8]*B.data[8];
		}
		public static void Prod(Mat3x3 A, Vec3 v, out Vec3 Ans) {
			Ans.x=A.data[0]*v.x+A.data[1]*v.y+A.data[2]*v.z;
			Ans.y=A.data[3]*v.x+A.data[4]*v.y+A.data[5]*v.z;
			Ans.z=A.data[6]*v.x+A.data[7]*v.y+A.data[8]*v.z;
		}
		public static void Prod(Vec3 v, Mat3x3 A, out Vec3 Ans) {
			Ans.x=A.data[0]*v.x+A.data[3]*v.y+A.data[6]*v.z;
			Ans.y=A.data[1]*v.x+A.data[4]*v.y+A.data[7]*v.z;
			Ans.z=A.data[2]*v.x+A.data[5]*v.y+A.data[8]*v.z;
		}
		public static void Invert(Mat3x3 A, Mat3x3 Ans) {
			double A11=A.data[4]*A.data[8]-A.data[5]*A.data[7], A21=A.data[5]*A.data[6]-A.data[3]*A.data[8], A31=A.data[3]*A.data[7]-A.data[4]*A.data[6];
			double det=A.data[0]*A11+A.data[1]*A21+A.data[2]*A31;
			Ans.data[0]=A11/det;
			Ans.data[1]=(A.data[7]*A.data[2]-A.data[1]*A.data[8])/det;
			Ans.data[2]=(A.data[1]*A.data[5]-A.data[4]*A.data[2])/det;
			Ans.data[3]=A21/det;
			Ans.data[4]=(A.data[0]*A.data[8]-A.data[6]*A.data[2])/det;
			Ans.data[5]=(-A.data[0]*A.data[5]+A.data[3]*A.data[2])/det;
			Ans.data[6]=A31/det;
			Ans.data[7]=(-A.data[0]*A.data[7]+A.data[6]*A.data[1])/det;
			Ans.data[8]=(A.data[0]*A.data[4]-A.data[3]*A.data[1])/det;
		}
		public static Mat3x3[] CreateArr(int leng) { Mat3x3[] ans=new Mat3x3[leng]; for (int i=0; i<leng; i++)ans[i]=new Mat3x3(); return ans; }
	}
	public class StdMat3x3 {
		static public void SetRotationMatrix(Mat3x3 A, Quaternion q) {
			double norm2=q.NormSqr();
			norm2=(norm2>0.0)?(2/norm2):0;
			double nx=q.x*norm2, ny=q.y*norm2, nz=q.z*norm2;
			double wx=q.w*nx, wy=q.w*ny, wz=q.w*nz;
			double xx=q.x*nx, yy=q.y*ny, zz=q.z*nz;
			double xy=q.x*ny, xz=q.x*nz, yz=q.y*nz;
			A.data[0]=1-yy-zz; A.data[1]=xy-wz; A.data[2]=xz+wy;
			A.data[3]=xy+wz; A.data[4]=1-xx-zz; A.data[5]=yz-wx;
			A.data[6]=xz-wy; A.data[7]=yz+wx; A.data[8]=1-xx-yy;
		}
		static public void SetRotationMatrixZXZ/*loc->abs Euler*/(Mat3x3 A, double s1, double c1, double s2, double c2, double s3, double c3) {
			A.data[0]=c1*c3-s1*c2*s3; A.data[1]=-c1*s3-s1*c2*c3; A.data[2]=s1*s2;
			A.data[3]=s1*c3+c1*c2*s3; A.data[4]=-s1*s3+c1*c2*c3; A.data[5]=-c1*s2;
			A.data[6]=s2*s3; A.data[7]=s2*c3; A.data[8]=c2;
		}
		static public void SetRotationMatrixZXZ/*loc->abs Euler*/(Mat3x3 A, double psi, double theta, double phi) {
			double c1=Math.Cos(psi), s1=Math.Sin(psi);
			double c2=Math.Cos(theta), s2=Math.Sin(theta);
			double c3=Math.Cos(phi), s3=Math.Sin(phi);
			SetRotationMatrixZXZ(A, s1, c1, s2, c2, s3, c3);
		}
		static public void SetRotationMatrixXYZ/*loc->abs Tait-Bryan*/(Mat3x3 A, double s1, double c1, double s2, double c2, double s3, double c3) {
			A.data[0]=c2*c3; A.data[1]=-c2*s3; A.data[2]=s2;
			A.data[3]=s1*s2*c3+c1*s3; A.data[4]=-s1*s2*s3+c1*c3; A.data[5]=-s1*c2;
			A.data[6]=-c1*s2*c3+s1*s3; A.data[7]=c1*s2*s3+s1*c3; A.data[8]=c1*c2;
		}
		static public void SetRotationMatrixXYZ/*loc->abs Tait-Bryan*/(Mat3x3 A, double alpha, double beta, double gamma) {
			double c1=Math.Cos(alpha), s1=Math.Sin(alpha);
			double c2=Math.Cos(beta), s2=Math.Sin(beta);
			double c3=Math.Cos(gamma), s3=Math.Sin(gamma);
			SetRotationMatrixXYZ(A, s1, c1, s2, c2, s3, c3);
		}
		static public void SetHessianMatrixXYZ/*loc->abs Tait-Bryan*/(Mat3x3 A, Mat3x3 A_1, Mat3x3 A_2, Mat3x3 A_3, Mat3x3 A_11, Mat3x3 A_22, Mat3x3 A_33, Mat3x3 A_12, Mat3x3 A_13, Mat3x3 A_23, double s1, double c1, double s2, double c2, double s3, double c3) {
			double t1=s1*s2, t2=t1*c3, t3=t2+c1*s3, t4=t1*s3, t5=-t4+c1*c3, t6=c1*s2, t7=t6*c3;
			double t8=-t7+s1*s3, t9=t6*s3, t10=t9+s1*c3, t11=c1*c2, t12=s1*c2, t13=c2*s3, t14=c2*c3;
			double t17=t12*s3, t18=t11*c3, t19=t12*c3, t20=t11*s3, t21=s2*c3, t22=s2*s3;
			A.Set(t14, -t13, s2, t3, t5, -t12, t8, t10, t11);
			A_1.Set(0, 0, 0, -t8, -t10, -t11, t3, t5, -t12);
			A_2.Set(-t21, t22, c2, t19, -t17, t1, -t18, t20, -t6);
			A_3.Set(-t13, -t14, 0, t5, -t3, 0, t10, -t8, 0);
			A_11.Set(0, 0, 0, -t3, -t5, t12, -t8, -t10, -t11);
			A_22.Set(-t14, t13, -s2, -t2, t4, t12, t7, -t9, -t11);
			A_33.Set(-t14, t13, 0, -t3, -t5, 0, -t8, -t10, 0);
			A_12.Set(0, 0, 0, t18, -t20, t6, t19, -t17, t1);
			A_13.Set(0, 0, 0, -t10, t8, 0, t5, -t3, 0);
			A_23.Set(t22, t21, 0, -t17, -t19, 0, t20, t18, 0);
		}
		static public void SetHessianMatrixXYZ/*loc->abs Tait-Bryan*/(Mat3x3 A, Mat3x3 A_1, Mat3x3 A_2, Mat3x3 A_3, Mat3x3 A_11, Mat3x3 A_22, Mat3x3 A_33, Mat3x3 A_12, Mat3x3 A_13, Mat3x3 A_23, double alpha, double beta, double gamma) {
			double c1=Math.Cos(alpha), s1=Math.Sin(alpha);
			double c2=Math.Cos(beta), s2=Math.Sin(beta);
			double c3=Math.Cos(gamma), s3=Math.Sin(gamma);
			SetHessianMatrixXYZ(A, A_1, A_2, A_3, A_11, A_22, A_33, A_12, A_13, A_23, s1, c1, s2, c2, s3, c3);
		}
		static public void GetRotationQuaternion(Mat3x3 A, ref Quaternion q) {
			q.w=0.5*Math.Sqrt(1+A.data[0]+A.data[4]+A.data[8]);
			if (q.w>1e-10) {
				double tmp2=0.25/q.w;
				q.x=(A.data[7]-A.data[5])*tmp2;
				q.y=(A.data[2]-A.data[6])*tmp2;
				q.z=(A.data[3]-A.data[1])*tmp2;
			} else {
				q.x=Math.Sign(A.data[7]-A.data[5])*Math.Sqrt(1+A.data[0]-A.data[4]-A.data[8])*0.5;
				q.y=Math.Sign(A.data[2]-A.data[6])*Math.Sqrt(1-A.data[0]+A.data[4]-A.data[8])*0.5;
				q.z=Math.Sign(A.data[3]-A.data[1])*Math.Sqrt(1-A.data[0]-A.data[4]+A.data[8])*0.5;
			}
		}
	}
	[Serializable]
	public struct Vec3 {
		public double x, y, z;
		public Vec3(double nx, double ny, double nz) { x=nx; y=ny; z=nz; }
		public Vec3(double[] xyz) { x=xyz[0]; y=xyz[1]; z=xyz[2]; }
		public Vec3(Vec3 v) { x=v.x; y=v.y; z=v.z; }
		public double[] Array { get { return new double[] { x, y, z }; } set { x=value[0]; y=value[1]; z=value[2]; } }
		public double this[int index] {
			get { switch (index) { case 0: { return x; } case 1: { return y; } case 2: { return z; } default: { return 0; } } }
			set { switch (index) { case 0: { x=value; break; } case 1: { y=value; break; } case 2: { z=value; break; } } }
		}
		public void Mul(double k) { x*=k; y*=k; z*=k; }
		public void Add(Vec3 v) { x+=v.x; y+=v.y; z+=v.z; }
		public void Sub(Vec3 v) { x-=v.x; y-=v.y; z-=v.z; }
		public void AA(double k, Vec3 v) { x=k*v.x; y=k*v.y; z=k*v.z; }
		public void MA(double k, Vec3 v) { x+=k*v.x; y+=k*v.y; z+=k*v.z; }
		public void Set(double nx, double ny, double nz) { x=nx; y=ny; z=nz; }
		public void Set(Vec3 v) { x=v.x; y=v.y; z=v.z; }
		public void Set0() { x=y=z=0; }
		public static Vec3 operator+(Vec3 v1, Vec3 v2) { return (new Vec3(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z)); }
		public static Vec3 operator-(Vec3 v1, Vec3 v2) { return (new Vec3(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z)); }
		public static Vec3 operator+(Vec3 v) { return (new Vec3(+v.x, +v.y, +v.z)); }
		public static Vec3 operator-(Vec3 v) { return (new Vec3(-v.x, -v.y, -v.z)); }
		public static Vec3 operator*(Vec3 v, double s) { return (new Vec3(v.x*s, v.y*s, v.z*s)); }
		public static Vec3 operator/(Vec3 v, double s) { return (new Vec3(v.x/s, v.y/s, v.z/s)); }
		public static Vec3 operator*(double s, Vec3 v) { return v*s; }
		public void Mul(Mat3x3 A) {
			double nx=A.data[0]*x+A.data[1]*y+A.data[2]*z;
			double ny=A.data[3]*x+A.data[4]*y+A.data[5]*z;
			double nz=A.data[6]*x+A.data[7]*y+A.data[8]*z;
			x=nx; y=ny; z=nz;
		}
		public static Vec3 CrossProd(Vec3 v1, Vec3 v2) { return (new Vec3(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x)); }
		public void CrossProduct(Vec3 v) { double tmp_x=x, tmp_y=y; x=tmp_y*v.z-z*v.y; y=z*v.x-tmp_x*v.z; z=tmp_x*v.y-tmp_y*v.x; }
		public static double ScalarProd(Vec3 v1, Vec3 v2) { return (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z); }
		public double ScalarProd(Vec3 v) { return (x*v.x+y*v.y+z*v.z); }
		public static double MixedProduct(Vec3 v1, Vec3 v2, Vec3 v3) { return v1.x*(v2.y*v3.z-v3.y*v2.z)+v1.y*(v2.z*v3.x-v2.x*v3.z)+v1.z*(v2.x*v3.y-v2.y*v3.x); }
		public static double Norm(Vec3 v) { return Math.Sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }
		public static double NormSqr(Vec3 v) { return v.x*v.x+v.y*v.y+v.z*v.z; }
		public static double NormS(Vec3 v) { return ((v.x>=0)?v.x:-v.x)+((v.z>=0)?v.z:-v.z)+((v.z>=0)?v.z:-v.z); }
		public double Norm() { return Math.Sqrt(x*x+y*y+z*z); }
		public double NormSqr() { return x*x+y*y+z*z; }
		public double NormS() { return ((x>=0)?x:-x)+((z>=0)?z:-z)+((z>=0)?z:-z); }
		public void Normalize() { double n=Norm(); if (Math.Abs(n)>1e-20)Mul(1/n); else Set0(); }
		public static Vec3 Interpolate(Vec3 v1, Vec3 v2, double lambda) { return (new Vec3(v1.x*(1-lambda)+v2.x*lambda, v1.y*(1-lambda)+v2.y*lambda, v1.z*(1-lambda)+v2.z*lambda)); }
		public static double Dist(Vec3 v1, Vec3 v2) { return (Math.Sqrt((v1.x-v2.x)*(v1.x-v2.x)+(v1.y-v2.y)*(v1.y-v2.y)+(v1.z-v2.z)*(v1.z-v2.z))); }
		public static double DistSqr(Vec3 v1, Vec3 v2) { return ((v1.x-v2.x)*(v1.x-v2.x)+(v1.y-v2.y)*(v1.y-v2.y)+(v1.z-v2.z)*(v1.z-v2.z)); }
		public static double DistS(Vec3 v1, Vec3 v2) { return Math.Abs(v1.x-v2.x)+Math.Abs(v1.y-v2.y)+Math.Abs(v1.z-v2.z); }
		public double Dist(Vec3 v) { return (Math.Sqrt((x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z))); }
		public double DistSqr(Vec3 v) { return ((x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z)); }
		public double DistS(Vec3 v) { return Math.Abs(x-v.x)+Math.Abs(y-v.y)+Math.Abs(z-v.z); }
		public static double Angle(Vec3 v1, Vec3 v2) { return Math.Acos(v1.ScalarProd(v2)/v1.Norm()/v2.Norm()); }
		public double Angle(Vec3 other) { return Angle(this, other); }
		public override string ToString() { return string.Format("(x={0,20} y={1,20} z={2,20})", x, y, z); }
		public static readonly Vec3 origin=new Vec3(0, 0, 0);
		public static readonly Vec3 AxisX=new Vec3(1, 0, 0);
		public static readonly Vec3 AxisY=new Vec3(0, 1, 0);
		public static readonly Vec3 AxisZ=new Vec3(0, 0, 1);
		public void Rotate(Quaternion q) {
			double t1=q.x*y-q.y*x+q.w*z;
			double t2=q.z*x-q.x*z+q.w*y;
			double t3=q.y*z-q.z*y+q.w*x;
			x+=2*(q.y*t1-q.z*t2);
			y+=2*(q.z*t3-q.x*t1);
			z+=2*(q.x*t2-q.y*t3);
		}
	}
	[Serializable]
	public struct Quaternion {
		public double w, x, y, z;
		public Quaternion(double nw, double nx, double ny, double nz) { x=nx; y=ny; z=nz; w=nw; }
		public Quaternion(double[] xyzw) { x=xyzw[1]; y=xyzw[2]; z=xyzw[3]; w=xyzw[0]; }
		public Quaternion(Quaternion v) { x=v.x; y=v.y; z=v.z; w=v.w; }
		public Quaternion(double nw) { x=y=z=0; w=nw; }
		public Quaternion(Vec3 v) { x=v.x; y=v.y; z=v.z; w=0; }
		public double[] Array { get { return new double[] { w, x, y, z }; } set { x=value[1]; y=value[2]; z=value[3]; w=value[0]; } }
		public double this[int index] {
			get {
				switch (index) {
				case 0: { return w; }
				case 1: { return x; }
				case 2: { return y; }
				case 3: { return z; }
				default: { return 0; }
				}
			}
			set {
				switch (index) {
				case 0: { w=value; break; }
				case 1: { x=value; break; }
				case 2: { y=value; break; }
				case 3: { z=value; break; }
				}
			}
		}
		public void SetEulerAngles(double psi, double theta, double phi) {
			double psi_c=Math.Cos(psi*0.5), psi_s=Math.Sin(psi*0.5);
			double phi_c=Math.Cos(phi*0.5), phi_s=Math.Sin(phi*0.5);
			double theta_c=Math.Cos(theta*0.5), theta_s=Math.Sin(theta*0.5);
			w=psi_c*theta_c*phi_c+phi_s*theta_s*phi_s;
			x=psi_c*theta_c*phi_s-phi_s*theta_s*phi_c;
			y=psi_c*theta_s*phi_c+phi_s*theta_c*phi_s;
			z=psi_s*theta_c*phi_c-phi_c*theta_s*phi_s;
		}
		public void GetEulerAngles(out double psi, out double theta, out double phi) {
			psi=Math.Atan2(2*(w*z+x*y), 1-2*(y*y+z*z));
			theta=Math.Asin(2*(w*y-x*z));
			phi=Math.Atan2(2*(w*x+y*z), 1-2*(x*x+y*y));
		}
		public void Mul(double k) { x*=k; y*=k; z*=k; w*=k; }
		public void Add(Quaternion v) { x+=v.x; y+=v.y; z+=v.z; w+=v.w; }
		public void AA(double k, Quaternion v) { x=k*v.x; y=k*v.y; z=k*v.z; w=k*v.w; }
		public void MA(double k, Quaternion v) { x+=k*v.x; y+=k*v.y; z+=k*v.z; w+=k*v.w; }
		public void Set(double nw, double nx, double ny, double nz) { x=nx; y=ny; z=nz; w=nw; }
		public void Set(Quaternion v) { x=v.x; y=v.y; z=v.z; w=v.w; }
		public void Set(Vec3 v) { x=v.x; y=v.y; z=v.z; w=0; }
		public void Set(double nw) { x=y=z=0; w=nw; }
		public void Set0() { x=y=z=w=0; }
		public void HamiltonProd(Quaternion v) {
			double lw=w, lx=x, ly=y, lz=z;
			w=lw*v.w-lx*v.x-ly*v.y-lz*v.z;
			x=lw*v.x+lx*v.w+ly*v.z-lz*v.y;
			y=lw*v.y-lx*v.z+ly*v.w+lz*v.x;
			z=lw*v.z+lx*v.y-ly*v.x+lz*v.w;
		}
		public void Conjugate() { x*=-1; y*=-1; z*=-1; }
		public Vec3 Rotate(Vec3 v) {
			double t1=-z*v.z-y*v.y-x*v.x;
			double t2=-z*v.y+y*v.z+w*v.x;
			double t3=z*v.x-x*v.z+w*v.y;
			double t4=-y*v.x+x*v.y+w*v.z;
			return new Vec3(w*t2-x*t1+y*t4-z*t3, w*t3-x*t4-y*t1+z*t2, w*t4+x*t3-y*t2-z*t1);
		}
		public Vec3 RotateBack(Vec3 v) {
			double t1=z*v.z+y*v.y+x*v.x;
			double t2=z*v.y-y*v.z+w*v.x;
			double t3=-z*v.x+x*v.z+w*v.y;
			double t4=y*v.x-x*v.y+w*v.z;
			return new Vec3(w*t2+x*t1-y*t4+z*t3, w*t3+x*t4+y*t1-z*t2, w*t4-x*t3+y*t2+z*t1);
		}
		public Quaternion Rotate(Quaternion v) {
			double t1=-z*v.z-y*v.y-x*v.x;
			double t2=-z*v.y+y*v.z+w*v.x;
			double t3=z*v.x-x*v.z+w*v.y;
			double t4=-y*v.x+x*v.y+w*v.z;
			return new Quaternion(v.w, w*t2-x*t1+y*t4-z*t3, w*t3-x*t4-y*t1+z*t2, w*t4+x*t3-y*t2-z*t1);
		}
		public Quaternion RotateBack(Quaternion v) {
			double t1=z*v.z+y*v.y+x*v.x;
			double t2=z*v.y-y*v.z+w*v.x;
			double t3=-z*v.x+x*v.z+w*v.y;
			double t4=y*v.x-x*v.y+w*v.z;
			return new Quaternion(v.w, w*t2+x*t1-y*t4+z*t3, w*t3+x*t4+y*t1-z*t2, w*t4-x*t3+y*t2+z*t1);
		}
		public static Quaternion HamiltonProd(Quaternion v1, Quaternion v2) { return new Quaternion(v1.w*v2.w-v1.x*v2.x-v1.y*v2.y-v1.z*v2.z, v1.w*v2.x+v1.x*v2.w+v1.y*v2.z-v1.z*v2.y, v1.w*v2.y-v1.x*v2.z+v1.y*v2.w+v1.z*v2.x, v1.w*v2.z+v1.x*v2.y-v1.y*v2.x+v1.z*v2.w); }
		public static Quaternion operator+(Quaternion v1, Quaternion v2) { return (new Quaternion(v1.w+v2.w, v1.x+v2.x, v1.y+v2.y, v1.z+v2.z)); }
		public static Quaternion operator-(Quaternion v1, Quaternion v2) { return (new Quaternion(v1.w-v2.w, v1.x-v2.x, v1.y-v2.y, v1.z-v2.z)); }
		public static Quaternion operator+(Quaternion v) { return (new Quaternion(+v.w, +v.x, +v.y, +v.z)); }
		public static Quaternion operator-(Quaternion v) { return (new Quaternion(-v.w, -v.x, -v.y, -v.z)); }
		public static Quaternion operator*(Quaternion v, double s) { return (new Quaternion(v.w*s, v.x*s, v.y*s, v.z*s)); }
		public static Quaternion operator/(Quaternion v, double s) { return (new Quaternion(v.w/s, v.x/s, v.y/s, v.z/s)); }
		public static Quaternion operator*(double s, Quaternion v) { return v*s; }
		public static Quaternion Conjugate(Quaternion q) { return new Quaternion(q.w, -q.x, -q.y, -q.z); }
		public static void Rotate(Quaternion q, ref Vec3 v) {
			double t1=-q.z*v.z-q.y*v.y-q.x*v.x;
			double t2=-q.z*v.y+q.y*v.z+q.w*v.x;
			double t3=q.z*v.x-q.x*v.z+q.w*v.y;
			double t4=-q.y*v.x+q.x*v.y+q.w*v.z;
			v.x=q.w*t2-q.x*t1+q.y*t4-q.z*t3;
			v.y=q.w*t3-q.x*t4-q.y*t1+q.z*t2;
			v.z=q.w*t4+q.x*t3-q.y*t2-q.z*t1;
		}
		public static void RotateBack(Quaternion q, ref Vec3 v) {
			double t1=q.z*v.z+q.y*v.y+q.x*v.x;
			double t2=q.z*v.y-q.y*v.z+q.w*v.x;
			double t3=-q.z*v.x+q.x*v.z+q.w*v.y;
			double t4=q.y*v.x-q.x*v.y+q.w*v.z;
			v.x=q.w*t2+q.x*t1-q.y*t4+q.z*t3;
			v.y=q.w*t3+q.x*t4+q.y*t1-q.z*t2;
			v.z=q.w*t4-q.x*t3+q.y*t2+q.z*t1;
		}
		public static void Rotate(Quaternion q, ref Quaternion v) {
			double t1=-q.z*v.z-q.y*v.y-q.x*v.x;
			double t2=-q.z*v.y+q.y*v.z+q.w*v.x;
			double t3=q.z*v.x-q.x*v.z+q.w*v.y;
			double t4=-q.y*v.x+q.x*v.y+q.w*v.z;
			v.x=q.w*t2-q.x*t1+q.y*t4-q.z*t3;
			v.y=q.w*t3-q.x*t4-q.y*t1+q.z*t2;
			v.z=q.w*t4+q.x*t3-q.y*t2-q.z*t1;
		}
		public static void RotateBack(Quaternion q, ref Quaternion v) {
			double t1=q.z*v.z+q.y*v.y+q.x*v.x;
			double t2=q.z*v.y-q.y*v.z+q.w*v.x;
			double t3=-q.z*v.x+q.x*v.z+q.w*v.y;
			double t4=q.y*v.x-q.x*v.y+q.w*v.z;
			v.x=q.w*t2+q.x*t1-q.y*t4+q.z*t3;
			v.y=q.w*t3+q.x*t4+q.y*t1-q.z*t2;
			v.z=q.w*t4-q.x*t3+q.y*t2+q.z*t1;
		}
		public static double ScalarProd(Quaternion v1, Quaternion v2) { return (v1.w*v2.w+v1.x*v2.x+v1.y*v2.y+v1.z*v2.z); }
		public double ScalarProd(Quaternion v) { return (w*v.w+x*v.x+y*v.y+z*v.z); }
		public double Dist(Quaternion v) { return (Math.Sqrt((w-v.w)*(w-v.w)+(x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z))); }
		public double DistSqr(Quaternion v) { return ((w-v.w)*(w-v.w)+(x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z)); }
		public double DistS(Quaternion v) { return Math.Abs(w-v.w)+Math.Abs(x-v.x)+Math.Abs(y-v.y)+Math.Abs(z-v.z); }
		public static double Dist(Quaternion v1, Quaternion v2) { return (Math.Sqrt((v1.w-v2.w)*(v1.w-v2.w)+(v1.x-v2.x)*(v1.x-v2.x)+(v1.y-v2.y)*(v1.y-v2.y)+(v1.z-v2.z)*(v1.z-v2.z))); }
		public static double DistSqr(Quaternion v1, Quaternion v2) { return ((v1.w-v2.w)*(v1.w-v2.w)+(v1.x-v2.x)*(v1.x-v2.x)+(v1.y-v2.y)*(v1.y-v2.y)+(v1.z-v2.z)*(v1.z-v2.z)); }
		public static double DistS(Quaternion v1, Quaternion v2) { return Math.Abs(v1.w-v2.w)+Math.Abs(v1.x-v2.x)+Math.Abs(v1.y-v2.y)+Math.Abs(v1.z-v2.z); }
		public static double Norm(Quaternion v) { return Math.Sqrt(v.w*v.w+v.x*v.x+v.y*v.y+v.z*v.z); }
		public static double NormSqr(Quaternion v) { return v.w*v.w+v.x*v.x+v.y*v.y+v.z*v.z; }
		public static double NormS(Quaternion v) { return ((v.w>=0)?v.w:-v.w)+((v.x>=0)?v.x:-v.x)+((v.z>=0)?v.z:-v.z)+((v.z>=0)?v.z:-v.z); }
		public double Norm() { return Math.Sqrt(w*w+x*x+y*y+z*z); }
		public double NormSqr() { return w*w+x*x+y*y+z*z; }
		public double NormS() { return ((w>=0)?w:-w)+((x>=0)?x:-x)+((z>=0)?z:-z)+((z>=0)?z:-z); }
		public void Normalize() { double n=Norm(); if (Math.Abs(n)>1e-20)Mul(1/n); else Set0(); }
		public static Quaternion Interpolate(Quaternion v1, Quaternion v2, double lambda) { return (new Quaternion(v1.w*(1-lambda)+v2.w*lambda, v1.x*(1-lambda)+v2.x*lambda, v1.y*(1-lambda)+v2.y*lambda, v1.z*(1-lambda)+v2.z*lambda)); }
		public override string ToString() { return string.Format("( w={0} x={1}, y={2}, z={3} ,size={4})", w, x, y, z, Norm()); }
		public static readonly Quaternion Id=new Quaternion(1, 0, 0, 0);
	}
	public class CoordSys {
		public Mat3x3 Transform_loc_to_abs=new Mat3x3();
		public Vec3 Pos_abs=Vec3.origin;
		public Vec3 Vel_abs=Vec3.origin;
		public Vec3 Vel_loc=Vec3.origin;
		public Quaternion Rotation=Quaternion.Id;
		public Vec3 Omega_abs=Vec3.origin;
		public Vec3 Omega_loc=Vec3.origin;
		public Vec3 AbsAxisX=Vec3.AxisX;
		public Vec3 AbsAxisY=Vec3.AxisY;
		public Vec3 AbsAxisZ=Vec3.AxisZ;
		public void ConvertVector_loc_to_abs(Vec3 loc_v, out Vec3 abs_v) { Mat3x3.Prod(Transform_loc_to_abs, loc_v, out abs_v); }
		public void ConvertVector_abs_to_loc(Vec3 abs_v, out Vec3 loc_v) { Mat3x3.Prod(abs_v, Transform_loc_to_abs, out loc_v); }
		public void ConvertPoint_loc_to_abs(Vec3 loc_p, out Vec3 abs_p) {
			Mat3x3.Prod(Transform_loc_to_abs, loc_p, out abs_p);
			abs_p.Add(Pos_abs);
		}
		public void ConvertPoint_abs_to_loc(Vec3 abs_p, out Vec3 loc_p) {
			Mat3x3.Prod(new Vec3(abs_p.x-Pos_abs.x, abs_p.y-Pos_abs.y, abs_p.z-Pos_abs.z), Transform_loc_to_abs, out loc_p);
		}
		public CoordSys() { }
		public CoordSys(CoordSys CS) {
			Transform_loc_to_abs=new Mat3x3(CS.Transform_loc_to_abs);
			Pos_abs=CS.Pos_abs;
			Vel_abs=CS.Vel_abs;
			Omega_abs=CS.Omega_abs;
			Vel_loc=CS.Vel_loc;
			Omega_loc=CS.Omega_loc;
			Rotation=CS.Rotation;
			AbsAxisX=CS.AbsAxisX;
			AbsAxisY=CS.AbsAxisY;
			AbsAxisZ=CS.AbsAxisZ;
		}
		public CoordSys(Vec3 pos_abs, Vec3 vel_abs, Quaternion rotation, Vec3 omega_loc) { Set(pos_abs, vel_abs, rotation, omega_loc); }
		public void Set(Vec3 pos_abs, Vec3 vel_abs, Quaternion rotation, Vec3 omega_loc) {
			Pos_abs=pos_abs;
			Vel_abs=vel_abs;
			Rotation=rotation/rotation.Norm();
			Omega_loc=omega_loc;
			StdMat3x3.SetRotationMatrix(Transform_loc_to_abs, Rotation);
			ConvertVector_abs_to_loc(Vel_abs, out Vel_loc);
			ConvertVector_loc_to_abs(Omega_loc, out Omega_abs);
			Transform_loc_to_abs.GetAxisVectors(out AbsAxisX, out AbsAxisY, out AbsAxisZ);
		}
		public void SetOmegaAbs(Vec3 pos_abs, Vec3 vel_abs, Quaternion rotation, Vec3 omega_abs) {
			Pos_abs=pos_abs;
			Vel_abs=vel_abs;
			Rotation=rotation/rotation.Norm();
			Omega_abs=omega_abs;
			StdMat3x3.SetRotationMatrix(Transform_loc_to_abs, Rotation);
			ConvertVector_abs_to_loc(Vel_abs, out Vel_loc);
			ConvertVector_abs_to_loc(Omega_abs, out Omega_loc);
			Transform_loc_to_abs.GetAxisVectors(out AbsAxisX, out AbsAxisY, out AbsAxisZ);
		}
		public void Get(out Vec3 pos_abs, out Vec3 vel_abs, out Quaternion rotation, out Vec3 omega_loc) {
			pos_abs=Pos_abs;
			vel_abs=Vel_abs;
			rotation=Rotation;
			omega_loc=Omega_loc;
		}
		public void GetFirstDerivations(ref Vec3 D_pos_abs, ref Quaternion D_rot) {
			D_pos_abs.Set(Vel_abs);
			D_rot.Set(
			0.5*(-Rotation.x*Omega_loc.x-Rotation.y*Omega_loc.y-Rotation.z*Omega_loc.z),
			0.5*(+Rotation.w*Omega_loc.x+Rotation.y*Omega_loc.z-Rotation.z*Omega_loc.y),
			0.5*(+Rotation.w*Omega_loc.y-Rotation.x*Omega_loc.z+Rotation.z*Omega_loc.x),
			0.5*(+Rotation.w*Omega_loc.z+Rotation.x*Omega_loc.y-Rotation.y*Omega_loc.x));
		}
	}
}
//namespace Utils {
//	using System;
//	public class Matrix3x3 {
//		public double[] data=new double[9];
//		public double this[int x, int y] { get { return data[x+3*y]; } set { data[x+3*y]=value; } }
//		public Matrix3x3() { Set0(); }
//		public Matrix3x3(Matrix3x3 m) { for(int i=0; i<9; i++)data[i]=m.data[i]; }
//		public void Mul(double k) { for(int i=0; i<9; i++)data[i]*=k; }
//		public void Add(Matrix3x3 m) { for(int i=0; i<9; i++)data[i]+=m.data[i]; }
//		public void AA(double k, Matrix3x3 m) { for(int i=0; i<9; i++)data[i]=k*m.data[i]; }
//		public void MA(double k, Matrix3x3 m) { for(int i=0; i<9; i++)data[i]+=k*m.data[i]; }
//		public void Set(Matrix3x3 m) { for(int i=0; i<9; i++)data[i]=m.data[i]; }
//		public void Set0() { for(int i=0; i<9; i++)data[i]=0; }
//		public void GetAxisVectors(out Vector3 Vx, out Vector3 Vy, out Vector3 Vz) { Vx.x=data[0]; Vy.x=data[1]; Vz.x=data[2]; Vx.y=data[3]; Vy.y=data[4]; Vz.y=data[5]; Vx.z=data[6]; Vy.z=data[7]; Vz.z=data[8]; }
//		public static void Prod(Matrix3x3 A, Matrix3x3 B, Matrix3x3 Ans) {
//			Ans.data[0]=A.data[0]*B.data[0]+A.data[1]*B.data[3]+A.data[2]*B.data[6];
//			Ans.data[1]=A.data[0]*B.data[1]+A.data[1]*B.data[4]+A.data[2]*B.data[7];
//			Ans.data[2]=A.data[0]*B.data[2]+A.data[1]*B.data[5]+A.data[2]*B.data[8];
//			Ans.data[3]=A.data[3]*B.data[0]+A.data[4]*B.data[3]+A.data[5]*B.data[6];
//			Ans.data[4]=A.data[3]*B.data[1]+A.data[4]*B.data[4]+A.data[5]*B.data[7];
//			Ans.data[5]=A.data[3]*B.data[2]+A.data[4]*B.data[5]+A.data[5]*B.data[8];
//			Ans.data[6]=A.data[6]*B.data[0]+A.data[7]*B.data[3]+A.data[8]*B.data[6];
//			Ans.data[7]=A.data[6]*B.data[1]+A.data[7]*B.data[4]+A.data[8]*B.data[7];
//			Ans.data[8]=A.data[6]*B.data[2]+A.data[7]*B.data[5]+A.data[8]*B.data[8];
//		}
//		public static void Prod(Matrix3x3 A, Vector3 v, out Vector3 Ans) {
//			Ans.x=A.data[0]*v.x+A.data[1]*v.y+A.data[2]*v.z;
//			Ans.y=A.data[3]*v.x+A.data[4]*v.y+A.data[5]*v.z;
//			Ans.z=A.data[6]*v.x+A.data[7]*v.y+A.data[8]*v.z;
//		}
//		public static void Prod(Vector3 v, Matrix3x3 A, out Vector3 Ans) {
//			Ans.x=A.data[0]*v.x+A.data[3]*v.y+A.data[6]*v.z;
//			Ans.y=A.data[1]*v.x+A.data[4]*v.y+A.data[7]*v.z;
//			Ans.z=A.data[2]*v.x+A.data[5]*v.y+A.data[8]*v.z;
//		}
//		public void SetRotationMatrix(Quaternion q) {
//			double norm2=q.NormSqr();
//			norm2=(norm2>0.0)?(2/norm2):0;
//			double nx=q.x*norm2, ny=q.y*norm2, nz=q.z*norm2;
//			double wx=q.w*nx, wy=q.w*ny, wz=q.w*nz;
//			double xx=q.x*nx, yy=q.y*ny, zz=q.z*nz;
//			double xy=q.x*ny, xz=q.x*nz, yz=q.y*nz;
//			data[0]=1-yy-zz; data[1]=xy-wz; data[2]=xz+wy;
//			data[3]=xy+wz; data[4]=1-xx-zz; data[5]=yz-wx;
//			data[6]=xz-wy; data[7]=yz+wx; data[8]=1-xx-yy;
//		}
//		public void SetRotationMatrix/*loc->abs*/(double psi, double theta, double phi) {
//			double psi_c=Math.Cos(psi*0.5), psi_s=Math.Sin(psi*0.5);
//			double phi_c=Math.Cos(phi*0.5), phi_s=Math.Sin(phi*0.5);
//			double theta_c=Math.Cos(theta*0.5), theta_s=Math.Sin(theta*0.5);
//			data[0]=psi_c*phi_c-psi_s*theta_c*phi_s;
//			data[1]=-psi_c*phi_s-psi_s*theta_c*phi_c;
//			data[2]=psi_s*theta_s;
//			data[3]=psi_s*phi_c+psi_c*theta_c*phi_s;
//			data[4]=-psi_s*phi_s+psi_c*theta_c*phi_c;
//			data[5]=-psi_c*theta_s;
//			data[6]=theta_s*phi_s;
//			data[7]=theta_s*phi_c;
//			data[8]=theta_c;
//		}
//		public void Transpose() {
//			double tmp=data[1]; data[1]=data[3]; data[3]=tmp;
//			tmp=data[2]; data[2]=data[6]; data[6]=tmp;
//			tmp=data[5]; data[5]=data[7]; data[7]=tmp;
//		}
//		public void GetRotationQuaternion(Quaternion q) {
//			double tmp1=Math.Sqrt(1+data[0]+data[4]+data[8]);
//			q.w=0.5*tmp1;
//			if(q.w>1e-10) {
//				double tmp2=0.5/tmp1;
//				q.x=(data[7]-data[5])*tmp2;
//				q.y=(data[2]-data[6])*tmp2;
//				q.z=(data[3]-data[1])*tmp2;
//			} else {
//				q.x=Math.Sign(data[7]-data[5])*Math.Sqrt(1+data[0]-data[4]-data[8])*0.5;
//				q.y=Math.Sign(data[2]-data[6])*Math.Sqrt(1-data[0]+data[4]-data[8])*0.5;
//				q.z=Math.Sign(data[3]-data[1])*Math.Sqrt(1-data[0]-data[4]+data[8])*0.5;
//			}
//		}
//	}
//	public struct Vector3 {
//		public double x, y, z;
//		public Vector3(double nx, double ny, double nz) { x=nx; y=ny; z=nz; }
//		public Vector3(double[] xyz) { x=xyz[0]; y=xyz[1]; z=xyz[2]; }
//		public Vector3(Vector3 v) { x=v.x; y=v.y; z=v.z; }
//		public double[] Array { get { return new double[] { x, y, z }; } set { x=value[0]; y=value[1]; z=value[2]; } }
//		public double this[int index] {
//			get {
//				switch(index) {
//					case 0: { return x; }
//					case 1: { return y; }
//					case 2: { return z; }
//					default: { return 0; }
//				}
//			}
//			set {
//				switch(index) {
//					case 0: { x=value; break; }
//					case 1: { y=value; break; }
//					case 2: { z=value; break; }
//				}
//			}
//		}
//		public void Mul(double k) { x*=k; y*=k; z*=k; }
//		public void Add(Vector3 v) { x+=v.x; y+=v.y; z+=v.z; }
//		public void AA(double k, Vector3 v) { x=k*v.x; y=k*v.y; z=k*v.z; }
//		public void MA(double k, Vector3 v) { x+=k*v.x; y+=k*v.y; z+=k*v.z; }
//		public void Set(double nx, double ny, double nz) { x=nx; y=ny; z=nz; }
//		public void Set(Vector3 v) { x=v.x; y=v.y; z=v.z; }
//		public void Set0() { x=y=z=0; }
//		public static Vector3 operator+(Vector3 v1, Vector3 v2) { return (new Vector3(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z)); }
//		public static Vector3 operator-(Vector3 v1, Vector3 v2) { return (new Vector3(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z)); }
//		public static Vector3 operator+(Vector3 v) { return (new Vector3(+v.x, +v.y, +v.z)); }
//		public static Vector3 operator-(Vector3 v) { return (new Vector3(-v.x, -v.y, -v.z)); }
//		public static Vector3 operator*(Vector3 v, double s) { return (new Vector3(v.x*s, v.y*s, v.z*s)); }
//		public static Vector3 operator/(Vector3 v, double s) { return (new Vector3(v.x/s, v.y/s, v.z/s)); }
//		public static Vector3 operator*(double s, Vector3 v) { return v*s; }
//		public void Prod(Matrix3x3 A) {
//			double nx=A.data[0]*x+A.data[1]*y+A.data[2]*z;
//			double ny=A.data[3]*x+A.data[4]*y+A.data[5]*z;
//			double nz=A.data[6]*x+A.data[7]*y+A.data[8]*z;
//			x=nx; y=ny; z=nz;
//		}
//		public static Vector3 CrossProd(Vector3 v1, Vector3 v2) { return (new Vector3(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x)); }
//		public void CrossProduct(Vector3 v) { double tmp_x=x, tmp_y=y; x=tmp_y*v.z-z*v.y; y=z*v.x-tmp_x*v.z; z=tmp_x*v.y-tmp_y*v.x; }
//		public static double ScalarProd(Vector3 v1, Vector3 v2) { return (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z); }
//		public double ScalarProd(Vector3 v) { return (x*v.x+y*v.y+z*v.z); }
//		public static double MixedProduct(Vector3 v1, Vector3 v2, Vector3 v3) { return v1.x*(v2.y*v3.z-v3.y*v2.z)+v1.y*(v2.z*v3.x-v2.x*v3.z)+v1.z*(v2.x*v3.y-v2.y*v3.x); }
//		public static double Norm(Vector3 v) { return Math.Sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }
//		public static double NormSqr(Vector3 v) { return v.x*v.x+v.y*v.y+v.z*v.z; }
//		public static double NormS(Vector3 v) { return ((v.x>=0)?v.x:-v.x)+((v.z>=0)?v.z:-v.z)+((v.z>=0)?v.z:-v.z); }
//		public double Norm() { return Math.Sqrt(x*x+y*y+z*z); }
//		public double NormSqr() { return x*x+y*y+z*z; }
//		public double NormS() { return ((x>=0)?x:-x)+((z>=0)?z:-z)+((z>=0)?z:-z); }
//		public void Normalize() { double n=Norm(); if(Math.Abs(n)>1e-20)Mul(1/n); else Set0(); }
//		public static Vector3 Interpolate(Vector3 v1, Vector3 v2, double lambda) { return (new Vector3(v1.x*(1-lambda)+v2.x*lambda, v1.y*(1-lambda)+v2.y*lambda, v1.z*(1-lambda)+v2.z*lambda)); }
//		public static double Dist(Vector3 v1, Vector3 v2) { return (Math.Sqrt((v1.x-v2.x)*(v1.x-v2.x)+(v1.y-v2.y)*(v1.y-v2.y)+(v1.z-v2.z)*(v1.z-v2.z))); }
//		public static double DistSqr(Vector3 v1, Vector3 v2) { return ((v1.x-v2.x)*(v1.x-v2.x)+(v1.y-v2.y)*(v1.y-v2.y)+(v1.z-v2.z)*(v1.z-v2.z)); }
//		public static double DistS(Vector3 v1, Vector3 v2) { return Math.Abs(v1.x-v2.x)+Math.Abs(v1.y-v2.y)+Math.Abs(v1.z-v2.z); }
//		public double Dist(Vector3 v) { return (Math.Sqrt((x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z))); }
//		public double DistSqr(Vector3 v) { return ((x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z)); }
//		public double DistS(Vector3 v) { return Math.Abs(x-v.x)+Math.Abs(y-v.y)+Math.Abs(z-v.z); }
//		public static double Angle(Vector3 v1, Vector3 v2) { return Math.Acos(v1.ScalarProd(v2)/v1.Norm()/v2.Norm()); }
//		public double Angle(Vector3 other) { return Angle(this, other); }
//		public override string ToString() { return string.Format("(x={0,20} y={1,20} z={2,20})", x, y, z); }
//		public static readonly Vector3 origin=new Vector3(0, 0, 0);
//		public static readonly Vector3 AxisX=new Vector3(1, 0, 0);
//		public static readonly Vector3 AxisY=new Vector3(0, 1, 0);
//		public static readonly Vector3 AxisZ=new Vector3(0, 0, 1);
//		public void Rotate(Quaternion q) {
//			double t1=q.x*y-q.y*x+q.w*z;
//			double t2=q.z*x-q.x*z+q.w*y;
//			double t3=q.y*z-q.z*y+q.w*x;
//			x+=2*(q.y*t1-q.z*t2);
//			y+=2*(q.z*t3-q.x*t1);
//			z+=2*(q.x*t2-q.y*t3);
//		}
//	}
//	public struct Quaternion {
//		public double w, x, y, z;
//		public Quaternion(double nw, double nx, double ny, double nz) { x=nx; y=ny; z=nz; w=nw; }
//		public Quaternion(double[] xyzw) { x=xyzw[1]; y=xyzw[2]; z=xyzw[3]; w=xyzw[0]; }
//		public Quaternion(Quaternion v) { x=v.x; y=v.y; z=v.z; w=v.w; }
//		public Quaternion(double nw) { x=y=z=0; w=nw; }
//		public Quaternion(Vector3 v) { x=v.x; y=v.y; z=v.z; w=0; }
//		public double[] Array { get { return new double[] { w, x, y, z }; } set { x=value[1]; y=value[2]; z=value[3]; w=value[0]; } }
//		public double this[int index] {
//			get {
//				switch(index) {
//					case 0: { return w; }
//					case 1: { return x; }
//					case 2: { return y; }
//					case 3: { return z; }
//					default: { return 0; }
//				}
//			}
//			set {
//				switch(index) {
//					case 0: { w=value; break; }
//					case 1: { x=value; break; }
//					case 2: { y=value; break; }
//					case 3: { z=value; break; }
//				}
//			}
//		}
//		public void SetEulerAngles(double psi, double theta, double phi) {
//			double psi_c=Math.Cos(psi*0.5), psi_s=Math.Sin(psi*0.5);
//			double phi_c=Math.Cos(phi*0.5), phi_s=Math.Sin(phi*0.5);
//			double theta_c=Math.Cos(theta*0.5), theta_s=Math.Sin(theta*0.5);
//			w=psi_c*theta_c*phi_c+phi_s*theta_s*phi_s;
//			x=psi_c*theta_c*phi_s-phi_s*theta_s*phi_c;
//			y=psi_c*theta_s*phi_c+phi_s*theta_c*phi_s;
//			z=psi_s*theta_c*phi_c-phi_c*theta_s*phi_s;
//		}
//		public void GetEulerAngles(out double psi, out double theta, out double phi) {
//			psi=Math.Atan2(2*(w*z+x*y), 1-2*(y*y+z*z));
//			theta=Math.Asin(2*(w*y-x*z));
//			phi=Math.Atan2(2*(w*x+y*z), 1-2*(x*x+y*y));
//		}
//		public void Mul(double k) { x*=k; y*=k; z*=k; w*=k; }
//		public void Add(Quaternion v) { x+=v.x; y+=v.y; z+=v.z; w+=v.w; }
//		public void AA(double k, Quaternion v) { x=k*v.x; y=k*v.y; z=k*v.z; w=k*v.w; }
//		public void MA(double k, Quaternion v) { x+=k*v.x; y+=k*v.y; z+=k*v.z; w+=k*v.w; }
//		public void Set(double nw, double nx, double ny, double nz) { x=nx; y=ny; z=nz; w=nw; }
//		public void Set(Quaternion v) { x=v.x; y=v.y; z=v.z; w=v.w; }
//		public void Set(Vector3 v) { x=v.x; y=v.y; z=v.z; w=0; }
//		public void Set(double nw) { x=y=z=0; w=nw; }
//		public void Set0() { x=y=z=w=0; }
//		public void HamiltonProd(Quaternion v) {
//			double lw=w, lx=x, ly=y, lz=z;
//			w=lw*v.w-lx*v.x-ly*v.y-lz*v.z;
//			x=lw*v.x+lx*v.w+ly*v.z-lz*v.y;
//			y=lw*v.y-lx*v.z+ly*v.w+lz*v.x;
//			z=lw*v.z+lx*v.y-ly*v.x+lz*v.w;
//		}
//		public static Quaternion HamiltonProd(Quaternion v1, Quaternion v2) { return new Quaternion(v1.w*v2.w-v1.x*v2.x-v1.y*v2.y-v1.z*v2.z, v1.w*v2.x+v1.x*v2.w+v1.y*v2.z-v1.z*v2.y, v1.w*v2.y-v1.x*v2.z+v1.y*v2.w+v1.z*v2.x, v1.w*v2.z+v1.x*v2.y-v1.y*v2.x+v1.z*v2.w); }
//		public static Quaternion operator+(Quaternion v1, Quaternion v2) { return (new Quaternion(v1.w+v2.w, v1.x+v2.x, v1.y+v2.y, v1.z+v2.z)); }
//		public static Quaternion operator-(Quaternion v1, Quaternion v2) { return (new Quaternion(v1.w-v2.w, v1.x-v2.x, v1.y-v2.y, v1.z-v2.z)); }
//		public static Quaternion operator+(Quaternion v) { return (new Quaternion(+v.w, +v.x, +v.y, +v.z)); }
//		public static Quaternion operator-(Quaternion v) { return (new Quaternion(-v.w, -v.x, -v.y, -v.z)); }
//		public static Quaternion operator*(Quaternion v, double s) { return (new Quaternion(v.w*s, v.x*s, v.y*s, v.z*s)); }
//		public static Quaternion operator/(Quaternion v, double s) { return (new Quaternion(v.w/s, v.x/s, v.y/s, v.z/s)); }
//		public static Quaternion operator*(double s, Quaternion v) { return v*s; }
//		public static double ScalarProd(Quaternion v1, Quaternion v2) { return (v1.w*v2.w+v1.x*v2.x+v1.y*v2.y+v1.z*v2.z); }
//		public double ScalarProd(Quaternion v) { return (w*v.w+x*v.x+y*v.y+z*v.z); }
//		public static double Norm(Quaternion v) { return Math.Sqrt(v.w*v.w+v.x*v.x+v.y*v.y+v.z*v.z); }
//		public static double NormSqr(Quaternion v) { return v.w*v.w+v.x*v.x+v.y*v.y+v.z*v.z; }
//		public static double NormS(Quaternion v) { return ((v.w>=0)?v.w:-v.w)+((v.x>=0)?v.x:-v.x)+((v.z>=0)?v.z:-v.z)+((v.z>=0)?v.z:-v.z); }
//		public double Norm() { return Math.Sqrt(w*w+x*x+y*y+z*z); }
//		public double NormSqr() { return w*w+x*x+y*y+z*z; }
//		public double NormS() { return ((w>=0)?w:-w)+((x>=0)?x:-x)+((z>=0)?z:-z)+((z>=0)?z:-z); }
//		public void Normalize() { double n=Norm(); if(Math.Abs(n)>1e-20)Mul(1/n); else Set0(); }
//		public static Quaternion Interpolate(Quaternion v1, Quaternion v2, double lambda) { return (new Quaternion(v1.w*(1-lambda)+v2.w*lambda, v1.x*(1-lambda)+v2.x*lambda, v1.y*(1-lambda)+v2.y*lambda, v1.z*(1-lambda)+v2.z*lambda)); }
//		public override string ToString() { return string.Format("( w={0} x={1}, y={2}, z={3} ,size={4})", w, x, y, z, Norm()); }
//		public static readonly Quaternion Id=new Quaternion(1, 0, 0, 0);
//	}
//	public class CoordSys {
//		public Matrix3x3 Transform_loc_to_abs=new Matrix3x3();
//		public Vector3 Pos_abs=Vector3.origin;
//		public Vector3 Vel_abs=Vector3.origin;
//		public Vector3 Vel_loc=Vector3.origin;
//		public Quaternion Rotation=Quaternion.Id;
//		public Vector3 Omega_abs=Vector3.origin;
//		public Vector3 Omega_loc=Vector3.origin;
//		public Vector3 AbsAxisX=Vector3.AxisX;
//		public Vector3 AbsAxisY=Vector3.AxisY;
//		public Vector3 AbsAxisZ=Vector3.AxisZ;
//		public void ConvertVector_loc_to_abs(Vector3 loc_v, out Vector3 abs_v) { Matrix3x3.Prod(Transform_loc_to_abs, loc_v, out abs_v); }
//		public void ConvertVector_abs_to_loc(Vector3 abs_v, out Vector3 loc_v) { Matrix3x3.Prod(abs_v, Transform_loc_to_abs, out loc_v); }
//		public void ConvertPoint_loc_to_abs(Vector3 loc_p, out Vector3 abs_p) {
//			Matrix3x3.Prod(Transform_loc_to_abs, loc_p, out abs_p);
//			abs_p.Add(Pos_abs);
//		}
//		public void ConvertPoint_abs_to_loc(Vector3 abs_p, out Vector3 loc_p) {
//			Matrix3x3.Prod(new Vector3(abs_p.x-Pos_abs.x, abs_p.y-Pos_abs.y, abs_p.z-Pos_abs.z), Transform_loc_to_abs, out loc_p);
//		}
//		public CoordSys() { }
//		public CoordSys(CoordSys CS) {
//			Transform_loc_to_abs=new Matrix3x3(CS.Transform_loc_to_abs);
//			Pos_abs=CS.Pos_abs;
//			Vel_abs=CS.Vel_abs;
//			Omega_abs=CS.Omega_abs;
//			Vel_loc=CS.Vel_loc;
//			Omega_loc=CS.Omega_loc;
//			Rotation=CS.Rotation;
//			AbsAxisX=CS.AbsAxisX;
//			AbsAxisY=CS.AbsAxisY;
//			AbsAxisZ=CS.AbsAxisZ;
//		}
//		public CoordSys(Vector3 pos_abs, Vector3 vel_abs, Quaternion rotation, Vector3 omega_loc) { Set(pos_abs, vel_abs, rotation, omega_loc); }
//		public void Set(Vector3 pos_abs, Vector3 vel_abs, Quaternion rotation, Vector3 omega_loc) {
//			Pos_abs=pos_abs;
//			Vel_abs=vel_abs;
//			Rotation=rotation/rotation.Norm();
//			Omega_loc=omega_loc;
//			Transform_loc_to_abs.SetRotationMatrix(Rotation);
//			ConvertVector_abs_to_loc(Vel_abs, out Vel_loc);
//			ConvertVector_loc_to_abs(Omega_loc, out Omega_abs);
//			Transform_loc_to_abs.GetAxisVectors(out AbsAxisX, out AbsAxisY, out AbsAxisZ);
//		}
//		public void SetOmegaAbs(Vector3 pos_abs, Vector3 vel_abs, Quaternion rotation, Vector3 omega_abs) {
//			Pos_abs=pos_abs;
//			Vel_abs=vel_abs;
//			Rotation=rotation/rotation.Norm();
//			Omega_abs=omega_abs;
//			Transform_loc_to_abs.SetRotationMatrix(Rotation);
//			ConvertVector_abs_to_loc(Vel_abs, out Vel_loc);
//			ConvertVector_abs_to_loc(Omega_abs, out Omega_loc);
//			Transform_loc_to_abs.GetAxisVectors(out AbsAxisX, out AbsAxisY, out AbsAxisZ);
//		}
//		public void Get(out Vector3 pos_abs, out Vector3 vel_abs, out Quaternion rotation, out Vector3 omega_loc) {
//			pos_abs=Pos_abs;
//			vel_abs=Vel_abs;
//			rotation=Rotation;
//			omega_loc=Omega_loc;
//		}
//		public void GetFirstDerivations(ref Vector3 D_pos_abs, ref Quaternion D_rot) {
//			D_pos_abs.Set(Vel_abs);
//			D_rot.Set(
//			0.5*(-Rotation.x*Omega_loc.x-Rotation.y*Omega_loc.y-Rotation.z*Omega_loc.z),
//			0.5*(+Rotation.w*Omega_loc.x+Rotation.y*Omega_loc.z-Rotation.z*Omega_loc.y),
//			0.5*(+Rotation.w*Omega_loc.y-Rotation.x*Omega_loc.z+Rotation.z*Omega_loc.x),
//			0.5*(+Rotation.w*Omega_loc.z+Rotation.x*Omega_loc.y-Rotation.y*Omega_loc.x));
//		}
//	}
//	class WWW : Car.World {
//		static public Vector3 GG=new Vector3(0, 0, 0);
//		static public Vector3 NN=new Vector3(0, 0, 1);//new Vector3(0.1736481777, 0, 0.9848077530);
//		private Vector3 tmp=Vector3.origin;
//		public bool RayCast(Vector3 point, Vector3 dir, double max_dist, out double dist, out double K_friction, out Vector3 Contact_Point, out Vector3 Normal, out string MaterialName) {
//			K_friction=1;
//			Normal=NN;
//			tmp.Set(GG); tmp.MA(-1, point);
//			double SurfaceDist=Vector3.ScalarProd(tmp, NN);
//			double DistCoeff=Vector3.ScalarProd(dir, NN);
//			dist=SurfaceDist/DistCoeff;
//			Contact_Point=point+dist*dir;
//			MaterialName="none";
//			return (dist>0)&&(dist<max_dist);
//		}
//		public WWW() { }
//	}
//}