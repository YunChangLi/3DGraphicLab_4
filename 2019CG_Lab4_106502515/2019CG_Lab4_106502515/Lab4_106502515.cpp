#include <gl/GLUT.H>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>

#define TMSize 4
#define PI 3.141592f

using namespace std;

//cd C:\Users\bruce\source\repos\2019CG_Lab4_106502515\Debug
//2019CG_Lab4_106502515.exe



struct ASC {

	int vertex;					//the number of the vertex
	int face;					//the number of the face
	float vertexList[4][4000];  //record the every point in the model
	float faceRGB[3][4000];		//record the every color in the face
	float faceEQ[4][4000];		//record the every plane equation in the face
	vector< vector<int> > faceList; //record the every face in the model

	float Or, Og, Ob, Kd, Ks, N ;  //N : specular exponent

};
struct lightSource {

	float	Ipr, Ipg, Ipb; // Ipr , Ipg , Ipb
	float	Ix, Iy, Iz; //lightsource

};
struct winBuffer {
	float z = 10000.0f;
	float R, G, B;
};
void displayFunc();
void loadFile();
void readASC(string);
void readFile(string);
void reset(float a[][TMSize]);
void clearScreen();
void printMatrix(float matrix[][TMSize]);
void multi_Matrix(float a[][TMSize], float b[][TMSize]); //a*b = b
void cross_Matrix(float a[], float b[], float c[]); //aXb = b
float dot_Matrix(float a[], float b[]); //a.b
float positiveDot(float a[], float b[]);
void getNormalVector(vector<int>, int , int , float[] );
void translating(float, float, float, float a[][TMSize]);
void scaling(float, float, float, float a[][TMSize]);
void rotating(float, float, float, float a[][TMSize]);
void EM(float, float, float, float, float, float, float);
void GRM(float, float, float, float, float, float);
void mirror();
void TiltM(float);
void UnitizeVector(float a[]);
void PM(float, float, float);
void noFaceback(bool c[], int num);  //isSeeTheFace , ascModel_ID
vector< vector<float> > clipping(vector< vector<float> >);
vector< vector<float> > planeOne(vector<float>, vector<float>, vector< vector<float> >);
vector< vector<float> > planeTwo(vector<float>, vector<float>, vector< vector<float> >);
vector< vector<float> > planeThr(vector<float>, vector<float>, vector< vector<float> >);
vector< vector<float> > planeFor(vector<float>, vector<float>, vector< vector<float> >);
vector< vector<float> > planeFiv(vector<float>, vector<float>, vector< vector<float> >);
vector< vector<float> > planeSix(vector<float>, vector<float>, vector< vector<float> >);
void drawFace(float a[][8000], int , float R , float G , float B );  //the vertices of face , the number of vertices 
void drawLine(int, int, float a[][8000], float R, float G, float B); // record the ID of two vertices , and model ID
void lineAlgorithm(int x0, int y0, int x1, int y1, float R, float G, float B);
void setpixel(int x, int y , float R , float G , float B);
void Color(vector<int> , int faceID , int modelID);
void planeEQ(vector<int>, int faceID, int modelID , float V_List[][8000]);
float getColor(char RGB , int lightNum , int modelID , float[] , float[]);
void getLightVector(int num_light, float origin[], float Vector[]);
void getReflectVector(float NVector[], float LVector[], float RVector[]);
void getViewVector(float origin[], float Vector[]);
void setColorInWin(float a[][8000] , int , int modelID , int faceID , bool isForward[4000] , float b[4000] , float N[], float D);
void setBgColorInWin(float V_List[][4]);
bool isLeftPoint(int , int , int , int , int , int);
bool isRightPoint(int, int, int, int, int, int);


char *file;
ifstream inputFile;
ifstream inputASC;
int ScreenHeight, ScreenWidth;

ASC ascModel[10];
lightSource lighting[10];

int num_ascModel = 0;
float vxl, vxr, vyb, vyt;
float AspectRatio;
float fatt;
float KaIar, KaIag, KaIab;
float Br, Bg, Bb;  //background RGB
int index;
float Ir = 0.0f, Ig = 1.0f, Ib = 1.0f;

float eyePos[3];
float transform_Matrix[TMSize][TMSize];
float eye_Matrix[TMSize][TMSize];
float perspective_Matrix[TMSize][TMSize];
float WTM[TMSize][TMSize];
winBuffer Buffer[600][600]; //store the every color in the window


int main(int argc, char *argv[]) {
	system("pause");

	file = argv[1];
	reset(transform_Matrix);
	reset(eye_Matrix);
	reset(perspective_Matrix);
	reset(WTM);

	inputFile.open(file);
	if (inputFile.fail()) {
		cout << "fail" << endl;
	}
	else {
		inputFile >> ScreenWidth;
		inputFile >> ScreenHeight;
	}
	inputFile.close();
	//cout << "size : " << ScreenWidth << " " << ScreenHeight << endl;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(ScreenWidth, ScreenHeight);
	glutCreateWindow("Your First GLUT Window!");
	glutDisplayFunc(displayFunc);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, ScreenWidth, 0.0, ScreenHeight);



	glFlush();
	glutMainLoop();



	return 0;
}
void displayFunc() {
	loadFile();
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(1.0);
	glFlush();
}

void loadFile() {
	string temp;
	inputFile.open(file);

	if (inputFile.fail()) {
		cout << "fail" << endl;
	}
	else {


		while (!inputFile.eof()) {
			inputFile >> temp;
			//cout << temp << endl;
			readFile(temp);
			temp = "";
			//system("pause");
		}
		inputFile.close();

	}

}
void readFile(string code) {
	string str;
	string ascName;
	float Sx, Sy, Sz;
	float Tx, Ty, Tz;
	float Ex, Ey, Ez;
	float COIx, COIy, COIz;
	float Tilt;
	float Hither, Yon, Hav;
	float Rx, Ry, Rz;

	//Index : multiple light source ; Ix,y,z : the light source(L)

	//R(reflect vector) : 2|(Normal좩L)|*Normal-L
	//L(light vector) : drawing point to light position
	//V(view vector) : drawing point to eye position 
	//the formular for lighting : 
	/*
	Ir = (KaIar)*Or + [(i = 1 to num )Kd*Ipr(i)*(Nor좩L(i))*Or] + [(i = 1 to num )Ks*Ipr(i)*(R(i)좩V)^N]
	Ig = (KaIag)*Og + [(i = 1 to num )Kd*Ipg(i)*(Nor좩L(i))*Og] + [(i = 1 to num )Ks*Ipg(i)*(R(i)좩V)^N]
	Ib = (KaIar)*Ob + [(i = 1 to num )Kd*Ipb(i)*(Nor좩L(i))*Ob] + [(i = 1 to num )Ks*Ipb(i)*(R(i)좩V)^N]
	*/


	if (code == "#") {
		getline(inputFile, str);
	}
	else if (code == "scale") {

		inputFile >> Sx;
		inputFile >> Sy;
		inputFile >> Sz;
		/*cout << "Sx : " << Sx << endl;
		cout << "Sy : " << Sy << endl;
		cout << "Sz : " << Sz << endl;
		cout << "scaling" << endl;*/

		scaling(Sx, Sy, Sz, transform_Matrix);
		//cout << "transform_Matrix" << endl;
		//printMatrix(transform_Matrix);
	}
	else if (code == "rotate") {
		inputFile >> Rx;
		inputFile >> Ry;
		inputFile >> Rz;
		/*cout << "Rx : " << Rx << endl;
		cout << "Ry : " << Ry << endl;
		cout << "Rz : " << Rz << endl;
		cout << "rotating" << endl;*/
		
		rotating(Rx, Ry, Rz, transform_Matrix);
		//cout << "transform_Matrix" << endl;
		//printMatrix(transform_Matrix);
	}
	else if (code == "translate") {
		inputFile >> Tx;
		inputFile >> Ty;
		inputFile >> Tz;
		/*cout << "Tx : " << Tx << endl;
		cout << "Ty : " << Ty << endl;
		cout << "Tz : " << Tz << endl;
		cout << "translating" << endl;*/

		translating(Tx, Ty, Tz, transform_Matrix);
		//cout << "transform_Matrix" << endl;
		//printMatrix(transform_Matrix);
		
	}
	else if (code == "ambient") {

		inputFile >> KaIar; inputFile >> KaIag; inputFile >> KaIab;
		//cout << KaIar << " " << KaIag << " " << KaIab << endl;
	}
	else if (code == "object") {
		inputFile >> ascName;
		//cout << ascName << endl;
		
		readASC(ascName);

		inputFile >> ascModel[num_ascModel].Or; 
		inputFile >> ascModel[num_ascModel].Og;
		inputFile >> ascModel[num_ascModel].Ob;
		inputFile >> ascModel[num_ascModel].Kd;
		inputFile >> ascModel[num_ascModel].Ks;
		inputFile >> ascModel[num_ascModel].N;
		/*cout<< ascModel[num_ascModel].Or<<" "<< ascModel[num_ascModel].Og<<" "<< ascModel[num_ascModel].Ob<<endl;
		cout << ascModel[num_ascModel].Kd << " " << ascModel[num_ascModel].Ks << endl << "N : "<<ascModel[num_ascModel].N << endl;*/
		
		
		
		//do the translation for the model
		float sum = 0;
		vector<float> tempList;

		for (int j = 0; j < ascModel[num_ascModel].vertex; j++) {   //transformation
			for (int k = 0; k < TMSize; k++) {
				sum = 0;
				for (int m = 0; m < TMSize; m++) {
					sum += transform_Matrix[k][m] * ascModel[num_ascModel].vertexList[m][j];
				}
				tempList.push_back(sum);
			}
			ascModel[num_ascModel].vertexList[0][j] = tempList[0]; //x
			ascModel[num_ascModel].vertexList[1][j] = tempList[1]; //y
			ascModel[num_ascModel].vertexList[2][j] = tempList[2]; //z
			ascModel[num_ascModel].vertexList[3][j] = tempList[3]; //w
			tempList.clear();
			//cout << ascModel[num_ascModel].vertexList[0][j] << " , " << ascModel[num_ascModel].vertexList[1][j] << " , " << ascModel[num_ascModel].vertexList[2][j] << " , " << ascModel[num_ascModel].vertexList[3][j] << endl;
		}
		

		num_ascModel++;
		//cout << "modelnum : " << num_ascModel << endl;

	}
	else if (code == "observer") {
		inputFile >> Ex; inputFile >> Ey; inputFile >> Ez;  //10 , 10 , -10
		inputFile >> COIx; inputFile >> COIy; inputFile >> COIz; //0 , 0 , 0
		inputFile >> Tilt; // 1
		inputFile >> Hither; inputFile >> Yon; inputFile >> Hav;

		/*cout << Ex << " "; cout << Ey << " "; cout << Ez << endl;
		cout << COIx << " "; cout << COIy << " "; cout << COIz << endl;
		cout << Tilt << endl;;
		cout << Hither << " "; cout << Yon << " "; cout << Hav << endl;*/

		EM(Ex, Ey, Ez, COIx, COIy, COIz, Tilt);
		eyePos[0] = Ex; eyePos[1] = Ey; eyePos[2] = Ez;
		//cout << "eye_Matrix-------------------------------------------------------" << endl;
		//printMatrix(eye_Matrix);

		PM(Hither, Yon, Hav);
		//cout << "perspective_Matrix-----------------------------------------------" << endl;
		//printMatrix(perspective_Matrix);

	}
	else if (code == "background") {
		inputFile >> Br; inputFile >> Bg; inputFile >> Bb;
	}
	else if (code == "light") {
		inputFile >> index; 
		inputFile >> lighting[index - 1].Ipr; 
		inputFile >> lighting[index - 1].Ipg; 
		inputFile >> lighting[index - 1].Ipb;
		inputFile >> lighting[index - 1].Ix;
		inputFile >> lighting[index - 1].Iy;
		inputFile >> lighting[index - 1].Iz;

		/*cout << index << endl ;
		cout << lighting[index - 1].Ipr << " " << lighting[index - 1].Ipg << " " << lighting[index - 1].Ipb << endl;
		cout << lighting[index - 1].Ix << " " << lighting[index - 1].Iy << " " << lighting[index - 1].Iz << endl;*/

		
	}
	else if (code == "viewport") {

		inputFile >> vxl ; inputFile >> vxr ; inputFile >> vyb ; inputFile >> vyt ;
		/*cout << "vxl : " << vxl << endl;
		cout << "vxr : " << vxr << endl;
		cout << "vyb : " << vyb << endl;
		cout << "vyt : " << vyt << endl;*/

		float wxl, wxr, wyb, wyt;

		wxl = wyb = -1;
		wxr = wyt = 1;

		perspective_Matrix[1][1] = (vxr - vxl) / (vyt - vyb); //Aspect Ratio
		AspectRatio = (vxr - vxl) / (vyt - vyb);
		float scalarX = (vxr - vxl) / (wxr - wxl);
		float scalarY = (vyt - vyb) / (wyt - wyb);
		float moveX = (vxl - wxl) * ScreenWidth / (wxr - wxl);
		float moveY = (vyb - wyb) * ScreenHeight / (wyt - wyb);
		float RatioX = ScreenWidth * scalarX / (wxr - wxl);
		float RatioY = ScreenHeight * scalarY / (wyt - wyb);

		float tempWTM[TMSize][TMSize] = { {RatioX , 0 , 0  , RatioX + moveX} ,
										{0 , RatioY , 0 , RatioY + moveY} ,
										{0 , 0 , 1 , 0} ,
										{0 , 0 , 0 , 1} };


		for (int i = 0; i < TMSize; i++) {
			for (int j = 0; j < TMSize; j++) {
				WTM[j][i] = tempWTM[j][i];
			}
		}


		/*cout << "afterScalar" << endl;
		cout << "scree : " << RatioX << " , " << RatioY << endl;
		cout << "scalar : " << scalarX << " , " << scalarY << endl;
		cout << "move : " << moveX << " , " << moveY << endl;*/
		//printMatrix(WTM);

	}
	else if (code == "display") {

		clearScreen();   //first

		for (int i = 0; i < 600; i++) {
			for (int j = 0; j < 600; j++) {
				Buffer[i][j].z = 10000.0f;
				Buffer[i][j].R = 0.0f;
				Buffer[i][j].G = 0.0f;
				Buffer[i][j].B = 0.0f;
			}
		}



		float finalMulti_Matrix[TMSize][TMSize] = { {1 , 0 , 0 , 0} ,
													{0 , 1 , 0 , 0} ,
													{0 , 0 , 1 , 0} ,
													{0 , 0 , 0 , 1} };

		float bgPoint[TMSize][TMSize] = { {-1 , -1 , 1 , 1} ,
										  {-1 , 1 , -1 , 1} ,
										  {0 , 0 , 0 , 0} ,
										  {1 , 1 , 1 , 1} };

		float bgPointAfterWVM[TMSize][TMSize];
		float tempBGList[4];  //just the relay to build the BG matrix
		float sum = 0;

		//set BG Color
		//BackGround
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < TMSize; k++) {
				sum = 0;
				for (int m = 0; m < TMSize; m++) {
					sum += WTM[k][m] * bgPoint[m][j];
				}
				tempBGList[k] = sum;
			}
			bgPointAfterWVM[0][j] = tempBGList[0]; //x
			bgPointAfterWVM[1][j] = tempBGList[1]; //y
			bgPointAfterWVM[2][j] = tempBGList[2]; //z   
			bgPointAfterWVM[3][j] = tempBGList[3]; //w

		}
		setBgColorInWin(bgPointAfterWVM);


		multi_Matrix(eye_Matrix, finalMulti_Matrix);
		multi_Matrix(perspective_Matrix, finalMulti_Matrix);

		

		/*cout << "total model : " << num_ascModel << endl;
		cout << "total light : " << index << endl;*/

		for (int i = 0; i < num_ascModel; i++) {  //trans the model vertex to the screen space

			float temp_vertexList[4][4000];
			float temp_vertexClipping[4][8000];  //store the vertices that muliplied the WVM
			float tempList[4];  //just the relay to build the vertices matrix
			float HeadtempList[4000];
			bool isSeeTheFace[4000];


			memset(isSeeTheFace, true, 4000);    //initialize

			sum = 0;

			//no backFace------------------------------------------------------------------------------

			noFaceback(isSeeTheFace, i);

			//get the color----------------------------------------------------------------------------

			vector<int> faceIDarray;

			for (int d = 0; d < ascModel[i].face; d++) {

				for (int r = 0; r < ascModel[i].faceList[d].size(); r++) {

					faceIDarray.push_back( ascModel[i].faceList[d][r] - 1 );
					
					
				}
				//system("pause");
				Color(faceIDarray , d , i);
				faceIDarray.clear();

			}



			

			//multi EM PM + perspective divide-------------------------------------------------------------------------------
			for (int j = 0; j < ascModel[i].vertex; j++) {
				for (int k = 0; k < TMSize; k++) {
					sum = 0;
					for (int m = 0; m < TMSize; m++) {
						sum += finalMulti_Matrix[k][m] * ascModel[i].vertexList[m][j];
					}
					tempList[k] = sum;
				}
				temp_vertexList[0][j] = tempList[0];/// tempList[3]; //x
				temp_vertexList[1][j] = tempList[1];/// tempList[3]; //y
				temp_vertexList[2][j] = tempList[2];/// tempList[3]; //z   
				temp_vertexList[3][j] = tempList[3];/// tempList[3]; //w
			}


			//clipping---------------------------------------------------------------------------------
			vector< vector<float> > vec; //store the new point in the face
			vector<float> arrPoint;  //store the point ID temporarily

			for (int d = 0; d < ascModel[i].face; d++) {
				vec.clear();
				for (int r = 0; r < ascModel[i].faceList[d].size(); r++) {
					arrPoint.clear();
					int ID = ascModel[i].faceList[d][r] - 1;
					for (int h = 0; h < 4; h++) {
						arrPoint.push_back(temp_vertexList[h][ID]);
					}
					vec.push_back(arrPoint);
				}
				//clipping for every face

				//vec = clipping(vec);

				for (int e = 0; e < vec.size(); e++) {
					vec[e][0] = vec[e][0] / vec[e][3];
					vec[e][1] = vec[e][1] / vec[e][3];
					vec[e][2] = vec[e][2] / vec[e][3];
					vec[e][3] = vec[e][3] / vec[e][3];
					HeadtempList[d] = vec[e][2];
				}


				//multi WTM ------------------------------------------------------
				for (int j = 0; j < vec.size(); j++) {
					for (int k = 0; k < TMSize; k++) {
						sum = 0;
						for (int m = 0; m < TMSize; m++) {
							sum += WTM[k][m] * vec[j][m];
						}
						tempList[k] = sum;
					}
					temp_vertexClipping[0][j] = tempList[0]; //x
					temp_vertexClipping[1][j] = tempList[1]; //y
					temp_vertexClipping[2][j] = tempList[2]; //z   
					temp_vertexClipping[3][j] = tempList[3]; //w

				}
				
				//get normal vector

				float VerA[3], VerB[3];
				float normalVector[3];
				float dcoeff;

				for (int i = 0; i < 3; i++) {

					VerA[i] = temp_vertexClipping[i][2] - temp_vertexClipping[i][0];
					VerB[i] = temp_vertexClipping[i][1] - temp_vertexClipping[i][0];

				}
				cross_Matrix(VerA, VerB, normalVector);
				dcoeff = (normalVector[0] * temp_vertexClipping[0][0] + normalVector[1] * temp_vertexClipping[1][0] + normalVector[2] * temp_vertexClipping[2][0]);
				cout << normalVector[0] << " , " << normalVector[1] << " , " << normalVector[2] << " , "<<dcoeff<<endl;
				//draw Face for the Flat--------------------------------------------------------

				setColorInWin(temp_vertexClipping , vec.size(), i , d , isSeeTheFace , HeadtempList , normalVector , dcoeff);


				//draw the face OLD-------------------------------------------------------------
				//if (isSeeTheFace[d]) {

				//drawFace(temp_vertexClipping, vec.size() , ascModel[i].faceRGB[0][d] , ascModel[i].faceRGB[1][d] , ascModel[i].faceRGB[2][d]);
				

				//}
			}
		}

		

		//draw----------------------------------------------------------------------------------------

		for (int i = 0; i < ScreenHeight; i++) {
			for (int j = 0; j < ScreenWidth; j++) {
				setpixel(j, i, Buffer[i][j].R, Buffer[i][j].G, Buffer[i][j].B);
				
			}
		}
		//initialize
		for (int i = 0; i < 600; i++) {
			for (int j = 0; j < 600; j++) {
				Buffer[i][j].z = 10000.0f;
				Buffer[i][j].R = 0.0f;
				Buffer[i][j].G = 0.0f;
				Buffer[i][j].B = 0.0f;
			}
		}
		glFlush();
		system("pause");
	
	}
	else if (code == "reset") {
	    
		reset(transform_Matrix);
	}
	else if (code == "end") {
		exit(0);
	}

}
void readASC(string name) {
	inputASC.open(name);

	if (inputASC.fail()) {
		cout << "fail" << endl;
	}
	else {
		ASC model;
		inputASC >> model.vertex;
		inputASC >> model.face;
		//cout << model.vertex << " " << model.face << endl;
		
		//system("pause") ;
		int count;
		int pointID;
		for (int i = 0; i < model.vertex; i++) {

			inputASC >> model.vertexList[0][i]; //x
			inputASC >> model.vertexList[1][i]; //y
			inputASC >> model.vertexList[2][i]; //z
			model.vertexList[3][i] = 1;         //w
			//cout << setw(5) << model.vertexList[0][i] << setw(5) << model.vertexList[1][i] << setw(5) << model.vertexList[2][i] << setw(5) << model.vertexList[2][i] << endl;
		
		}
		
		for (int j = 0; j < model.face; j++) {
			inputASC >> count;
			//cout << count ;
			vector<int> v;
			for (int h = 0; h < count; h++) {
				inputASC >> pointID;
				//cout << setw(5) <<pointID ;
				v.push_back(pointID);
			}
			//cout << endl;
			model.faceList.push_back(v);
		}
		
		ascModel[num_ascModel] = model;
		inputASC.close();
		
	}
}
void printMatrix(float matrix[][TMSize]) {
	for (int i = 0; i < TMSize; i++) {
		cout << "[";
		for (int j = 0; j < TMSize; j++) {
			cout << setw(15) << matrix[i][j];
		}
		cout << "]" << endl;
	}
}
void multi_Matrix(float tm[][TMSize], float unit[][TMSize]) {

	float num = 0;
	float newMatrix[TMSize][TMSize];
	for (int i = 0; i < TMSize; i++) {
		for (int j = 0; j < TMSize; j++) {
			num = 0;
			for (int k = 0; k < TMSize; k++) {
				num += tm[i][k] * unit[k][j];
			}
			newMatrix[i][j] = num;
		}
	}

	for (int i = 0; i < TMSize; i++) {
		for (int j = 0; j < TMSize; j++) {
			unit[i][j] = newMatrix[i][j];
		}
	}

}
void cross_Matrix(float a[], float b[], float c[]) { //a X b return c
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
	//cout << "Cnormal " << a[1] * b[2] - a[2] * b[1] << " , " << a[2] * b[0] - a[0] * b[2] << " , " << a[0] * b[1] - a[1] * b[0] << endl;
}
float dot_Matrix(float a[], float b[]) {
	float result = 0;
	result += a[0] * b[0];
	result += a[1] * b[1];
	result += a[2] * b[2];

	return result;
}
float positiveDot(float a[], float b[]) {
	float result = 0;
	result += a[0] * b[0];
	result += a[1] * b[1];
	result += a[2] * b[2];
	if (result < 0) {
		result = 0;
	}

	return result;
}
void translating(float Tx, float Ty, float Tz, float a[][TMSize]) {
	float translate_Matrix[TMSize][TMSize] = { {1 , 0 , 0 , Tx},
												  {0 , 1 , 0 , Ty},
												  {0 , 0 , 1 , Tz},
												  {0 , 0 , 0 , 1} };

	multi_Matrix(translate_Matrix, a);
}
void scaling(float Sx, float Sy, float Sz, float a[][TMSize]) {
	float scale_Matrix[TMSize][TMSize] = { {Sx , 0 , 0 , 0},
											{0 , Sy , 0 , 0},
											{0 , 0 , Sz , 0},
											{0 , 0 , 0 , 1} };

	multi_Matrix(scale_Matrix, a);
}
void rotating(float Rx, float Ry, float Rz, float a[][TMSize]) {
	Rx = Rx * PI / 180.0;
	Ry = Ry * PI / 180.0;
	Rz = Rz * PI / 180.0;

	float rotateX_Matrix[TMSize][TMSize] = { {1 ,    0 ,         0      , 0},
											{0  , cos(Rx) , -1 * sin(Rx) , 0},
											{0  , sin(Rx) ,   cos(Rx)  , 0},
											{0 , 0 , 0 , 1} };

	float rotateY_Matrix[TMSize][TMSize] = { {cos(Ry) , 0 , sin(Ry) , 0},
											{0 ,       1 ,     0   , 0},
											{-1 * sin(Ry), 0 , cos(Ry) , 0},
											{0 , 0 , 0 , 1} };

	float rotateZ_Matrix[TMSize][TMSize] = { {cos(Rz) , -1 * sin(Rz) , 0 , 0},
											{sin(Rz) ,   cos(Rz) , 0 , 0},
											{0 , 0 , 1 ,0},
											{0 , 0 , 0 , 1} };

	multi_Matrix(rotateX_Matrix, a);
	multi_Matrix(rotateY_Matrix, a);
	multi_Matrix(rotateZ_Matrix, a);
}
void EM(float Ex, float Ey, float Ez, float COIx, float COIy, float COIz, float Tilt) {
	reset(eye_Matrix);
	translating(-1 * Ex, -1 * Ey, -1 * Ez, eye_Matrix);
	GRM(Ex, Ey, Ez, COIx, COIy, COIz);
	mirror();
	TiltM(Tilt);
}
void GRM(float Ex, float Ey, float Ez, float COIx, float COIy, float COIz) {
	//cout << "ExEyEz :  " << Ex << " " << Ey << " " << Ez << endl;
	float view_vector[3] = { COIx - Ex , COIy - Ey , COIz - Ez };
	float top_vector[3] = { 0 , 1 , 0 };

	float V1[3];
	float V2[3];
	float V3[3] = { COIx - Ex , COIy - Ey , COIz - Ez };


	UnitizeVector(V3);
	cross_Matrix(top_vector, view_vector, V1);   //get V1
	UnitizeVector(V1);
	cross_Matrix(V3, V1, V2);  //get V2
	UnitizeVector(V2);



	float GRM[TMSize][TMSize] = { {V1[0] , V1[1] , V1[2] , 0} ,
								 {V2[0] , V2[1] , V2[2] , 0} ,
								 {V3[0] , V3[1] , V3[2] , 0} ,
								 {0 , 0 , 0 , 1} };
	//cout << "GRM--------------------------------------------------------------------" << endl;
	//printMatrix(GRM);
	multi_Matrix(GRM, eye_Matrix);
}
void mirror() {
	float mirror_Matrix[TMSize][TMSize] = { {-1 , 0 , 0 , 0} ,
											{0 , 1 , 0 , 0} ,
											{0 , 0 , 1 , 0} ,
											{0 , 0 , 0 , 1} };

	//cout << "Mirror-----------------------------------------------------------------" << endl;
	//printMatrix(mirror_Matrix);
	multi_Matrix(mirror_Matrix, eye_Matrix);
}
void TiltM(float tilt) {
	tilt = tilt * PI / 180.0;
	float Tilt_Matrix[TMSize][TMSize] = { {cos(tilt)  , sin(tilt) , 0 , 0} ,
										{-1 * sin(tilt) , cos(tilt) , 0 , 0} ,
											{0 , 0 , 1 , 0} ,
											{0 , 0 , 0 , 1} };
	//cout << "Tilt--------------------------------------------------------------------" << endl;
	//printMatrix(Tilt_Matrix);
	multi_Matrix(Tilt_Matrix, eye_Matrix);
}
void PM(float H, float y, float theta) {

	theta = theta * PI / 180.0;
	float PM_22 = y * tan(theta) / (y - H);
	float PM_23 = H * y * tan(theta) / (H - y);
	float PM_32 = tan(theta);
	float PMtemp[TMSize][TMSize] = { {1 , 0 , 0 , 0} ,
									{0 , AspectRatio , 0 , 0} ,
									{0 , 0 , PM_22 , PM_23} ,
									{0 , 0 , PM_32 , 0} };


	for (int i = 0; i < TMSize; i++) {      //update the perspect Matrix
		for (int j = 0; j < TMSize; j++) {
			perspective_Matrix[i][j] = PMtemp[i][j];
		}
	}

}
void UnitizeVector(float vec[]) {
	float divisor = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	vec[0] /= divisor; vec[1] /= divisor; vec[2] /= divisor;
}
void noFaceback(bool isDrawFace[], int num) {

	for (int i = 0; i < ascModel[num].face; i++) {

		float Normal[3]; float EYE[3];

		if (ascModel[num].faceList[i].size() == 3) {  //triangle
			float VerA[3]; float VerB[3]; //VerA cross VerB		
			int ID_pointA = ascModel[num].faceList[i][0] - 1;  //call the point ID for cross
			int ID_pointB = ascModel[num].faceList[i][1] - 1;
			int ID_pointC = ascModel[num].faceList[i][2] - 1;

			for (int j = 0; j < 3; j++) {
				VerA[j] = ascModel[num].vertexList[j][ID_pointC] - ascModel[num].vertexList[j][ID_pointA];
				VerB[j] = ascModel[num].vertexList[j][ID_pointB] - ascModel[num].vertexList[j][ID_pointA];
				//calculate the EYE
				EYE[j] = ascModel[num].vertexList[j][ID_pointA] - eyePos[j];
			}
			cross_Matrix(VerA, VerB, Normal);  //get the normal vector
		}
		else if (ascModel[num].faceList[i].size() == 4) {   //square
			float VerA[3]; float VerB[3];
			int ID_pointA = ascModel[num].faceList[i][0] - 1;  //call the point ID for cross
			int ID_pointB = ascModel[num].faceList[i][1] - 1;
			int ID_pointC = ascModel[num].faceList[i][3] - 1;
			for (int j = 0; j < 3; j++) {
				VerA[j] = ascModel[num].vertexList[j][ID_pointC] - ascModel[num].vertexList[j][ID_pointA];
				VerB[j] = ascModel[num].vertexList[j][ID_pointB] - ascModel[num].vertexList[j][ID_pointA];
				//calculate the EYE
				EYE[j] = ascModel[num].vertexList[j][ID_pointA] - eyePos[j];
			}
			cross_Matrix(VerA, VerB, Normal);  //get the normal vector
		}
		if (dot_Matrix(Normal, EYE) > 0) {   //cannot see
			isDrawFace[i] = false;
		}


	}
}
vector< vector<float> > clipping(vector< vector<float> > sum) {
	//6 plane
	//a -- > b

	vector< vector<float> > fillsum;
	fillsum = sum;  //update

	vector< vector<float> > tempsum;
	for (int x = 0; x < fillsum.size(); x++) {  //plane One
		if (x == fillsum.size() - 1) {
			tempsum = planeOne(fillsum[x], fillsum[0], tempsum);
		}
		else {
			tempsum = planeOne(fillsum[x], fillsum[x + 1], tempsum);
		}
	}
	fillsum = tempsum;  //update
	tempsum.clear();

	for (int x = 0; x < fillsum.size(); x++) { //plane Sec
		if (x == fillsum.size() - 1) {
			tempsum = planeTwo(fillsum[x], fillsum[0], tempsum);
		}
		else {
			tempsum = planeTwo(fillsum[x], fillsum[x + 1], tempsum);
		}
	}
	fillsum = tempsum;  //update
	tempsum.clear();

	for (int x = 0; x < fillsum.size(); x++) {
		if (x == fillsum.size() - 1) {
			tempsum = planeThr(fillsum[x], fillsum[0], tempsum);
		}
		else {
			tempsum = planeThr(fillsum[x], fillsum[x + 1], tempsum);
		}
	}

	fillsum = tempsum;  //update
	tempsum.clear();

	for (int x = 0; x < fillsum.size(); x++) {
		if (x == fillsum.size() - 1) {
			tempsum = planeFor(fillsum[x], fillsum[0], tempsum);
		}
		else {
			tempsum = planeFor(fillsum[x], fillsum[x + 1], tempsum);
		}
	}

	fillsum = tempsum;  //update
	tempsum.clear();

	for (int x = 0; x < fillsum.size(); x++) {
		if (x == fillsum.size() - 1) {
			tempsum = planeFiv(fillsum[x], fillsum[0], tempsum);
		}
		else {
			tempsum = planeFiv(fillsum[x], fillsum[x + 1], tempsum);
		}
	}

	fillsum = tempsum;  //update
	tempsum.clear();

	for (int x = 0; x < fillsum.size(); x++) {
		if (x == fillsum.size() - 1) {
			tempsum = planeSix(fillsum[x], fillsum[0], tempsum);
		}
		else {
			tempsum = planeSix(fillsum[x], fillsum[x + 1], tempsum);
		}
	}

	fillsum = tempsum;  //update
	tempsum.clear();

	return fillsum;
}
vector< vector<float> > planeOne(vector<float> a, vector<float> b, vector< vector<float> > sum) { //clip w-x >= 0


	if (a[3] - a[0] >= 0 && b[3] - b[0] >= 0) {  // keep all
		sum.push_back(b);
	}
	else if (a[3] - a[0] < 0 && b[3] - b[0] < 0) { //reject all
		//nothing
	}
	else {

		//calculate the intersection
		float c1 = a[3] - a[0], c2 = b[3] - b[0];  //c = w-x
		float t = c1 / (c1 - c2);
		vector<float> inter;

		for (int i = 0; i < 4; i++) { //get the intersection point
			inter.push_back(a[i] + t * (b[i] - a[i]));
		}

		if (a[3] - a[0] < 0) {  //a : out ; b : in
			sum.push_back(inter);
			sum.push_back(b);

		}
		else if (b[3] - b[0] < 0) {  //a : in ; b : out
			sum.push_back(inter);
		}
	}
	return sum;
}
vector< vector<float> > planeTwo(vector<float> a, vector<float> b, vector< vector<float> > sum) { // clip W+X >= 0


	if (a[3] + a[0] >= 0 && b[3] + b[0] >= 0) {
		sum.push_back(b);
	}
	else if (a[3] + a[0] < 0 && b[3] + b[0] < 0) { //reject all

	}
	else {
		//calculate the intersection
		float c1 = a[3] + a[0], c2 = b[3] + b[0];  //c = w+x
		float t = c1 / (c1 - c2);

		vector<float> inter;
		for (int i = 0; i < 4; i++) { //get the intersection point
			inter.push_back(a[i] + t * (b[i] - a[i]));
		}

		if (a[3] + a[0] < 0) {  //a : out ; b : in
			sum.push_back(inter);
			sum.push_back(b);

		}
		else if (b[3] + b[0] < 0) {  //a : in ; b : out
			sum.push_back(inter);
		}

	}
	return sum;
}
vector< vector<float> > planeThr(vector<float> a, vector<float> b, vector< vector<float> > sum) { // clip W-Y >= 0


	if (a[3] - a[1] >= 0 && b[3] - b[1] >= 0) {
		sum.push_back(b);
	}
	else if (a[3] - a[1] < 0 && b[3] - b[1] < 0) { //reject all
		//nothing
	}
	else {
		//calculate the intersection
		float c1 = a[3] - a[1], c2 = b[3] - b[1];  //c = w-y
		float t = c1 / (c1 - c2);
		vector<float> inter;
		for (int i = 0; i < 4; i++) { //get the intersection point
			inter.push_back(a[i] + t * (b[i] - a[i]));
		}

		if (a[3] - a[1] < 0) {  //a : out ; b : in
			sum.push_back(inter);
			sum.push_back(b);

		}
		else if (b[3] - b[1] < 0) {  //a : in ; b : out
			sum.push_back(inter);
		}
	}
	return sum;
}
vector< vector<float> > planeFor(vector<float> a, vector<float> b, vector< vector<float> > sum) { // clip  W+Y >= 0


	if (a[3] + a[1] >= 0 && b[3] + b[1] >= 0) {
		sum.push_back(b);
	}
	else if (a[3] + a[1] < 0 && b[3] + b[1] < 0) { //reject all
		//nothing
	}
	else {
		//calculate the intersection
		float c1 = a[3] + a[1], c2 = b[3] + b[1];  //c = w+y
		float t = c1 / (c1 - c2);
		vector<float> inter;
		for (int i = 0; i < 4; i++) { //get the intersection point
			inter.push_back(a[i] + t * (b[i] - a[i]));
		}

		if (a[3] + a[1] < 0) {  //a : out ; b : in
			sum.push_back(inter);
			sum.push_back(b);

		}
		else if (b[3] + b[1] < 0) {  //a : in ; b : out
			sum.push_back(inter);
		}

	}
	return sum;
}
vector< vector<float> > planeFiv(vector<float> a, vector<float> b, vector< vector<float> > sum) { // clip w-z>=0


	if (a[3] - a[2] >= 0 && b[3] - b[2] >= 0) { //kepp all
		sum.push_back(b);
	}
	else if (a[3] - a[2] < 0 && b[3] - b[2] < 0) { //reject all
		//nothing
	}
	else {
		//calculate the intersection
		float c1 = a[3] - a[2], c2 = b[3] - b[2];  //c = w-z
		float t = c1 / (c1 - c2);
		vector<float> inter;
		for (int i = 0; i < 4; i++) { //get the intersection point
			inter.push_back(a[i] + t * (b[i] - a[i]));
		}

		if (a[3] - a[2] < 0) {  //a : out ; b : in
			sum.push_back(inter);
			sum.push_back(b);

		}
		else if (b[3] - b[2] < 0) {  //a : in ; b : out
			sum.push_back(inter);
		}

	}
	return sum;
}
vector< vector<float> > planeSix(vector<float> a, vector<float> b, vector< vector<float> > sum) { // clip z >= 0


	if (a[2] >= 0 && b[2] >= 0) {
		sum.push_back(b);
	}
	else if (a[2] < 0 && b[2] < 0) { //reject all
		//nothing
	}
	else {
		//calculate the intersection
		float c1 = a[2], c2 = b[2];  //c = z
		float t = c1 / (c1 - c2);
		vector<float> inter;
		for (int i = 0; i < 4; i++) { //get the intersection point
			inter.push_back(a[i] + t * (b[i] - a[i]));
		}

		if (a[2] < 0) {  //a : out ; b : in
			sum.push_back(inter);
			sum.push_back(b);

		}
		else if (b[2] < 0) {  //a : in ; b : out
			sum.push_back(inter);
		}
	}
	return sum;
}
void getNormalVector(vector<int> IDs,  int sizee, int modelID , float normal[]) {

	float VerA[3], VerB[3];
	
	
	for (int i = 0; i < 3; i++) {

		
		VerA[i] = ascModel[modelID].vertexList[i][IDs[2]] - ascModel[modelID].vertexList[i][IDs[0]];
		VerB[i] = ascModel[modelID].vertexList[i][IDs[1]] - ascModel[modelID].vertexList[i][IDs[0]];
		

	}
	/*cout <<"B : " <<VerB[0] << " , " << VerB[1] << " , " << VerB[2] << endl;
	cout << "A : " << VerA[0] << " , " << VerA[1] << " , " << VerA[2] << endl;*/
	cross_Matrix(VerA , VerB , normal);
}
void setpixel(int x, int y, float R, float G, float B)
{

	glBegin(GL_POINTS);/**start**/
	glColor3f(R ,G, B);  
	glVertex2i(x, y);
	glEnd();/**End**/
	return;

}
void drawFace(float V_List[][8000], int sizee , float R , float G , float B) {

	

	for (int i = 0; i < sizee; i++) {  //number of the vetices of a face

		if (i == sizee - 1) {  //the last one go back to the first one
			drawLine(i, 0, V_List , R , G , B);
		}
		else {
			drawLine(i, i + 1, V_List , R , G , B);
		}
	}

}
void drawLine(int cur, int next, float V_List[][8000], float R, float G, float B) {
	int x = V_List[0][cur];
	int y = V_List[1][cur];
	int xt = V_List[0][next];
	int yt = V_List[1][next];
	//cout << x << " , " << y << " , " << V_List[2][cur] << " ; " << xt << " , " << yt << " , " << V_List[2][next] << endl;
	//system("pause");
	lineAlgorithm(x, y, xt, yt , R , G , B);

}
void lineAlgorithm(int x0, int y0, int x1, int y1, float R , float G, float B)
{
	int firstx, firsty, secondx, secondy;
	firstx = x0; firsty = y0; secondx = x1; secondy = y1;
	/*cout << "x1 , y1 :" << firstx << " , " << firsty << endl;
	cout << "x2 , y2 :" << secondx << " , " << secondy << endl;*/

	/**initial judge**/
	bool inverse = abs(secondy - firsty) > abs(secondx - firstx);


	if (inverse) {   //swap(x1 , y1) swap(x2 , y2)  //scope > 1 || scope < -1
		int temp1 = 0;
		temp1 = firstx;
		firstx = firsty;
		firsty = temp1;

		temp1 = secondx;
		secondx = secondy;
		secondy = temp1;

	}
	if (firstx > secondx) //swap(x1 , x2) swap(y1 , y2)
	{
		int temp2 = 0;
		temp2 = firstx;
		firstx = secondx;
		secondx = temp2;

		temp2 = firsty;
		firsty = secondy;
		secondy = temp2;
	}


	/**initial**/
	bool Yminus = firsty > secondy;     //y1 < y2

	int x = firstx, y = firsty;
	int xt = secondx, yt = secondy;
	int a = yt - y, b = x - xt;
	int d = 0;
	if (Yminus) { d = a - b / 2; }
	else { d = a + b / 2; }


	while (x <= xt) {   //x <= x2

		if (inverse) {
			//if (overView(y, x)) {
			setpixel(y, x , R , G , B);
			//}
			//else { cout << "not print!!!" << endl; }

		}
		else {
			//if (overView(x, y)) {
			setpixel(x, y , R , G , B);
			//}
			//else { cout << "not print!!!" << endl; }
		}

		if (Yminus) {
			if (d >= 0) //choose E
			{
				x++;
				d += a;
			}
			else {      //SE(Yminus : 그쾤쨛⒀0)

				x++; y--;
				d += (a - b);  //SE(Yminus : 그쾤쨛⒀0)

			}
		}
		else {
			if (d <= 0) //choose E
			{

				x++;
				d += a;
			}
			else {      //choose NE 

				x++;  y++;
				d += (a + b);

			}
		}

	}

	return;
}
void Color(vector<int> IDs , int faceID , int modelID) {

	//get color for every face ---------------------------------------------------------

	//get the normal vector

	float normalVector[3];
	float originPoint[3];
	//get the normalVector and original point
	getNormalVector(IDs , IDs.size(), modelID , normalVector);
	UnitizeVector(normalVector);
	
	originPoint[0] = ascModel[modelID].vertexList[0][IDs[0]]; 
	originPoint[1] = ascModel[modelID].vertexList[1][IDs[0]];
	originPoint[2] = ascModel[modelID].vertexList[2][IDs[0]];


	Ir = getColor('r', index, modelID, normalVector, originPoint);
	Ig = getColor('g', index, modelID, normalVector, originPoint);
	Ib = getColor('b', index, modelID, normalVector, originPoint);

	//store Ir , Ig , Ib
	ascModel[modelID].faceRGB[0][faceID] = Ir;
	ascModel[modelID].faceRGB[1][faceID] = Ig;
	ascModel[modelID].faceRGB[2][faceID] = Ib;


}
void planeEQ(vector<int> IDs, int faceID, int modelID , float V_List[][8000]) {
	float normalVector[3];
	float originPoint[3];
	float d = 0.0f;
	
	//get the normalVector and original point
	//getNormalVector(IDs, IDs.size(), modelID, normalVector);
	float VerA[3], VerB[3];

	for (int i = 0; i < 3; i++) {

		VerA[i] = V_List[i][IDs[2]] - V_List[i][IDs[0]];
		VerB[i] = V_List[i][IDs[1]] - V_List[i][IDs[0]];

	}

	cross_Matrix(VerA, VerB, normalVector);

	//cout <<"Normal : "<<normalVector[0] << " , " << normalVector[1] << " , " << normalVector[2] << endl;
	originPoint[0] = V_List[0][IDs[0]];
	originPoint[1] = V_List[1][IDs[0]];
	originPoint[2] = V_List[2][IDs[0]];

	d = (normalVector[0] * originPoint[0] + normalVector[1] * originPoint[1] + normalVector[2] * originPoint[2]);
	

	//get the eq of plane
	ascModel[modelID].faceEQ[0][faceID] = normalVector[0];
	ascModel[modelID].faceEQ[1][faceID] = normalVector[1];
	ascModel[modelID].faceEQ[2][faceID] = normalVector[2];
	ascModel[modelID].faceEQ[3][faceID] = d ;

}
float getColor(char RGB , int lightNum , int modelID , float normal[] , float originPoint[]) {

	float reflectVector[3]; //R(reflect vector) : 2|(Normal좩L)| * Normal - L
	float lightVector[3];	//L(light vector) : drawing point to light position
	float viewVector[3];	//V(view vector) : drawing point to eye position

	float colorR = 0.0f, colorG = 0.0f, colorB = 0.0f;

	float diffuse = 0.0f;
	float specular = 0.0f;
	float ambient = 0.0f;

	getViewVector(originPoint, viewVector);
	UnitizeVector(viewVector);
	

	if (RGB == 'r') {
		//Ir = (KaIar)*Or + [(i = 1 to num)Kd*Ipr(i)*(Nor좩L(i))*Or] + [(i = 1 to num)Ks*Ipr(i)*(R(i)좩V) ^ N]
		
		ambient = KaIar * ascModel[modelID].Or;

		for (int i = 0; i < lightNum; i++) {

			getLightVector(i , originPoint ,lightVector ); //light number , origin , lightVector
			UnitizeVector(lightVector);

			getReflectVector(normal , lightVector ,reflectVector);  // check the light vector 
			UnitizeVector(reflectVector);
			

			diffuse += ascModel[modelID].Kd*lighting[i].Ipr * positiveDot(normal, lightVector)*ascModel[modelID].Or;
			specular += ascModel[modelID].Ks*lighting[i].Ipr * pow(positiveDot(reflectVector, viewVector), ascModel[modelID].N);
		}

		colorR = ambient + diffuse + specular;

		return colorR;
	}
	else if (RGB == 'g') {
		//Ig = (KaIag)*Og + [(i = 1 to num)Kd*Ipg(i)*(Nor좩L(i))*Og] + [(i = 1 to num)Ks*Ipg(i)*(R(i)좩V) ^ N]

		ambient = KaIag * ascModel[modelID].Og;

		for (int i = 0; i < lightNum; i++) {

			getLightVector(i, originPoint, lightVector); //light number , origin , lightVector
			UnitizeVector(lightVector);

			getReflectVector(normal, lightVector, reflectVector);
			UnitizeVector(reflectVector);
			

			diffuse += ascModel[modelID].Kd*lighting[i].Ipg * positiveDot(normal, lightVector)*ascModel[modelID].Og;
			specular += ascModel[modelID].Ks*lighting[i].Ipg * pow(positiveDot(reflectVector, viewVector), ascModel[modelID].N);
		}

		colorG = ambient + diffuse + specular;

		return colorG;
	}
	else if (RGB == 'b') {
		//Ib = (KaIar)*Ob + [(i = 1 to num)Kd*Ipb(i)*(Nor좩L(i))*Ob] + [(i = 1 to num)Ks*Ipb(i)*(R(i)좩V) ^ N]

		ambient = KaIab * ascModel[modelID].Ob;

		for (int i = 0; i < lightNum; i++) {

			getLightVector(i , originPoint ,lightVector ); //light number , origin , lightVector
			UnitizeVector(lightVector);

			getReflectVector(normal , lightVector ,reflectVector);
			UnitizeVector(reflectVector);
			

			diffuse += ascModel[modelID].Kd*lighting[i].Ipb * positiveDot(normal, lightVector)*ascModel[modelID].Ob;
			specular += ascModel[modelID].Ks*lighting[i].Ipb * pow(positiveDot(reflectVector, viewVector), ascModel[modelID].N);
		}

		colorB = ambient + diffuse + specular;

		return colorB;
	}
}
void getReflectVector(float NVector[], float LVector[], float RVector[]) {
	//2 | (Normal좩L) | *Normal - L
	float scalar;
	if (dot_Matrix(NVector, LVector) < 0) {
		scalar = 0;
	}
	else {
		scalar = 2 * dot_Matrix(NVector, LVector);
	}

	RVector[0] = scalar * NVector[0] - LVector[0];
	RVector[1] = scalar * NVector[1] - LVector[1];
	RVector[2] = scalar * NVector[2] - LVector[2];

}
void getLightVector(int num_light , float origin[], float Vector[]) {


	Vector[0] = lighting[num_light].Ix - origin[0];
	Vector[1] = lighting[num_light].Iy - origin[1];
	Vector[2] = lighting[num_light].Iz - origin[2];

}
void getViewVector(float origin[], float Vector[]) {
	Vector[0] = eyePos[0] - origin[0];
	Vector[1] = eyePos[1] - origin[1];
	Vector[2] = eyePos[2] - origin[2];
}
void setColorInWin(float V_List[][8000] , int sizee , int modelID , int faceID , bool isForward[4000] , float b[4000] , float N[] , float D) {

	/*get the xRange yRange size*/
	int minX = 10000, minY = 10000, maxX = -10000, maxY = -10000;

	for (int w = 0; w < sizee; w++) {
		
		if (V_List[0][w] < minX) { minX = V_List[0][w]; }
		if (V_List[0][w] > maxX) { maxX = V_List[0][w]; }
		if (V_List[1][w] < minY) { minY = V_List[1][w]; }
		if (V_List[1][w] > maxY) { maxY = V_List[1][w]; }
	}

	for (int i = minY; i <= maxY; i++) {
		for (int j = minX; j <= maxX; j++) {

			bool isInner = true;
			//cout << "Size : " << sizee << endl;
			for (int k = sizee - 1; k >= 0 ; k --) {

				int p1X, p1Y, p2X, p2Y;

				if (k == 0) {
					p1X = V_List[0][k]; p1Y = V_List[1][k];
					p2X = V_List[0][sizee - 1]; p2Y = V_List[1][sizee - 1];
				}
				else {
					p1X = V_List[0][k]; p1Y = V_List[1][k];
					p2X = V_List[0][k-1]; p2Y = V_List[1][k-1];
				}

				if (isForward[faceID]) {
					if (!isLeftPoint(p1Y, p1X, p2Y, p2X, i, j)) {

						isInner = false;
						break;
					}
				}
				else if(!isForward[faceID]){
					if (!isRightPoint(p1Y, p1X, p2Y, p2X, i, j)) {

						isInner = false;
						break;
					}
				}	

			}
			//system("pause");
			if (isInner) {

				float depthZ;
				depthZ = -1 * (N[0] * j + N[1] * i - D);

				depthZ /= N[2];

				
				if (depthZ < Buffer[i][j].z) {

					Buffer[i][j].R = ascModel[modelID].faceRGB[0][faceID];
					Buffer[i][j].G = ascModel[modelID].faceRGB[1][faceID];
					Buffer[i][j].B = ascModel[modelID].faceRGB[2][faceID];

					Buffer[i][j].z = depthZ;
				}

			}
		}
	}
}
void setBgColorInWin(float V_List[][4]) {
	/*get the xRange yRange size*/
	int minX = 10000, minY = 10000, maxX = -10000, maxY = -10000;

	for (int w = 0; w < 4; w++) {

		if (V_List[0][w] < minX) { minX = V_List[0][w]; }
		if (V_List[0][w] > maxX) { maxX = V_List[0][w]; }
		if (V_List[1][w] < minY) { minY = V_List[1][w]; }
		if (V_List[1][w] > maxY) { maxY = V_List[1][w]; }
	}

	for (int i = minY; i < maxY; i++) {
		for (int j = minX; j < maxX; j++) {


				float depthZ = 10000.0f;


				Buffer[i][j].R = Br;
				Buffer[i][j].G = Bg;
				Buffer[i][j].B = Bb;

				Buffer[i][j].z = depthZ;

		}
	}

}
bool isLeftPoint(int p1Y , int p1X , int p2Y , int p2X , int wY , int wX) {

	int tmpx = (p2Y - p1Y) * (wX - p1X) - (p2X - p1X) * (wY - p1Y);

	if (tmpx <= 0) {
		return true;
	}
	else { return false; }

}
bool isRightPoint(int p1Y, int p1X, int p2Y, int p2X, int wY, int wX) {

	int tmpx = (p2Y - p1Y) * (wX - p1X) - (p2X - p1X) * (wY - p1Y);

	if (tmpx >= 0) {
		return true;
	}
	else { return false; }

}
void reset(float matrix[][TMSize]) {
	for (int i = 0; i < TMSize; i++) {
		for (int j = 0; j < TMSize; j++) {
			if (i == j) {
				matrix[i][j] = 1;
			}
			else {
				matrix[i][j] = 0;
			}

		}
	}
}
void clearScreen() {

	glClear(GL_COLOR_BUFFER_BIT);
	//glFinish();
	glFlush();
}