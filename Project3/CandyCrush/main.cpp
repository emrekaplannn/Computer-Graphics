#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <GL/glew.h>   // The GL Header File
#include <GL/gl.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header



#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

GLuint gProgram;
int gWidth, gHeight;
int mposx=-1 , mposy=-1;
int objnumy=-1 , objnumx=-1;
bool clicked=false;


struct objeozellikleri
{
    objeozellikleri(bool a, bool b, bool c , float d, float e , int type_color){
        color = type_color;
        boom=a;
        buyuyor=b;
        kayiyor=c;
        scale=d;
        sliding=e;
    }
    int color;
    bool boom;
    bool buyuyor;
    bool kayiyor;
    float scale;
    float sliding;
};

std::vector<std::vector<objeozellikleri>> objeler;

struct Vertex
{
    Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Texture
{
    Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
    GLfloat u, v;
};

struct Normal
{
    Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Face
{
	Face(int v[], int t[], int n[]) {
		vIndex[0] = v[0];
		vIndex[1] = v[1];
		vIndex[2] = v[2];
		tIndex[0] = t[0];
		tIndex[1] = t[1];
		tIndex[2] = t[2];
		nIndex[0] = n[0];
		nIndex[1] = n[1];
		nIndex[2] = n[2];
	}
    GLuint vIndex[3], tIndex[3], nIndex[3];
};

vector<Vertex> gVertices;
vector<Texture> gTextures;
vector<Normal> gNormals;
vector<Face> gFaces;

GLuint gVertexAttribBuffer, gTextVBO, gIndexBuffer;
GLint gInVertexLoc, gInNormalLoc;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes;



bool ParseObj(const string& fileName)
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == '#') // comment
                {
                    continue;
                }
                else if (curLine[0] == 'v')
                {
                    if (curLine[1] == 't') // texture
                    {
                        str >> tmp; // consume "vt"
                        str >> c1 >> c2;
                        gTextures.push_back(Texture(c1, c2));
                    }
                    else if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gNormals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gVertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
					char c;
					int vIndex[3],  nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0]; 
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1]; 
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2]; 

					assert(vIndex[0] == nIndex[0] &&
						   vIndex[1] == nIndex[1] &&
						   vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

                    gFaces.push_back(Face(vIndex, tIndex, nIndex));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            //data += curLine;
            if (!myfile.eof())
            {
                //data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

	/*
	for (int i = 0; i < gVertices.size(); ++i)
	{
		Vector3 n;

		for (int j = 0; j < gFaces.size(); ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				if (gFaces[j].vIndex[k] == i)
				{
					// face j contains vertex i
					Vector3 a(gVertices[gFaces[j].vIndex[0]].x, 
							  gVertices[gFaces[j].vIndex[0]].y,
							  gVertices[gFaces[j].vIndex[0]].z);

					Vector3 b(gVertices[gFaces[j].vIndex[1]].x, 
							  gVertices[gFaces[j].vIndex[1]].y,
							  gVertices[gFaces[j].vIndex[1]].z);

					Vector3 c(gVertices[gFaces[j].vIndex[2]].x, 
							  gVertices[gFaces[j].vIndex[2]].y,
							  gVertices[gFaces[j].vIndex[2]].z);

					Vector3 ab = b - a;
					Vector3 ac = c - a;
					Vector3 normalFromThisFace = (ab.cross(ac)).getNormalized();
					n += normalFromThisFace;
				}

			}
		}

		n.normalize();

		gNormals.push_back(Normal(n.x, n.y, n.z));
	}
	*/

	assert(gVertices.size() == gNormals.size());

    return true;
}

bool ReadDataFromFile(
    const string& fileName, ///< [in]  Name of the shader file
    string&       data)     ///< [out] The contents of the file
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

void createVS()
{
    string shaderSource;

    string filename("vert.glsl");
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

    glAttachShader(gProgram, vs);
}

void createFS()
{
    string shaderSource;

    string filename("frag.glsl");
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

    glAttachShader(gProgram, fs);
}

void initShaders()
{
    gProgram = glCreateProgram();

    createVS();
    createFS();

    glLinkProgram(gProgram);
    glUseProgram(gProgram);
}

void initVBO()
{
    GLuint vao;
    glGenVertexArrays(1, &vao);
    assert(vao > 0);
    glBindVertexArray(vao);
    cout << "vao = " << vao << endl;

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	assert(glGetError() == GL_NONE);

	glGenBuffers(1, &gVertexAttribBuffer);
	glGenBuffers(1, &gIndexBuffer);

	assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
	gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
	int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);
	GLfloat* vertexData = new GLfloat [gVertices.size() * 3];
	GLfloat* normalData = new GLfloat [gNormals.size() * 3];
	GLuint* indexData = new GLuint [gFaces.size() * 3];

    float minX = 1e6, maxX = -1e6;
    float minY = 1e6, maxY = -1e6;
    float minZ = 1e6, maxZ = -1e6;

	for (int i = 0; i < gVertices.size(); ++i)
	{
		vertexData[3*i] = gVertices[i].x;
		vertexData[3*i+1] = gVertices[i].y;
		vertexData[3*i+2] = gVertices[i].z;

        minX = std::min(minX, gVertices[i].x);
        maxX = std::max(maxX, gVertices[i].x);
        minY = std::min(minY, gVertices[i].y);
        maxY = std::max(maxY, gVertices[i].y);
        minZ = std::min(minZ, gVertices[i].z);
        maxZ = std::max(maxZ, gVertices[i].z);
	}

    std::cout << "minX = " << minX << std::endl;
    std::cout << "maxX = " << maxX << std::endl;
    std::cout << "minY = " << minY << std::endl;
    std::cout << "maxY = " << maxY << std::endl;
    std::cout << "minZ = " << minZ << std::endl;
    std::cout << "maxZ = " << maxZ << std::endl;

	for (int i = 0; i < gNormals.size(); ++i)
	{
		normalData[3*i] = gNormals[i].x;
		normalData[3*i+1] = gNormals[i].y;
		normalData[3*i+2] = gNormals[i].z;
	}

	for (int i = 0; i < gFaces.size(); ++i)
	{
		indexData[3*i] = gFaces[i].vIndex[0];
		indexData[3*i+1] = gFaces[i].vIndex[1];
		indexData[3*i+2] = gFaces[i].vIndex[2];
	}


	glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

	// done copying; can free now
	delete[] vertexData;
	delete[] normalData;
	delete[] indexData;

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
}





void init(string s) 
{
	ParseObj(s);
	//ParseObj("bunny.obj");

    glEnable(GL_DEPTH_TEST);
    initShaders();
    initVBO();
}

void drawModel()
{
	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

	glDrawElements(GL_TRIANGLES, gFaces.size() * 3, GL_UNSIGNED_INT, 0);
}

void matchingFunction(int a, int b,int m){

    int patlayacak=0;
    //for rows;
    for(int i = 0 ; i<(a) ; i++){
        for(int y = 0  ; y<(b-2) ; y++){

            if(objeler[i][y].color == objeler[i][y+1].color && objeler[i][y+1].color == objeler[i][y+2].color){
                    objeler[i][y].buyuyor = true;
                    objeler[i][y+1].buyuyor = true;
                    objeler[i][y+2].buyuyor = true;
                    
            }
            if(objeler[i+1][y].color == objeler[i][y].color && objeler[i+1][y].color == objeler[i+2][y].color){
                    objeler[i+1][y].buyuyor = true;
                    objeler[i][y].buyuyor = true;
                    objeler[i+2][y].buyuyor = true;       
            }

        }
    }
    for(int i = 0 ; i<(b) ; i++){
        for(int y = 0  ; y<(a-2) ; y++){

            if(objeler[i+1][y].color == objeler[i][y].color && objeler[i+2][y].color == objeler[i][y].color){
                    objeler[i+1][y].buyuyor = true;
                    objeler[i][y].buyuyor = true;
                    objeler[i+2][y].buyuyor = true;
                    
            }

        }
    }
    
    

}


void display(int a, int b, float m,int toplam)
{   
    
    int counter = 0 ;
    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    
	static float angle = 0, slide_down=0, explode=1;
    
    int num_of_objx=a , num_of_objy=b;
    float mr= (float) 10/num_of_objx;
    float mr2= (float) 9/num_of_objy;
    
   
        for(int x=1; x<=num_of_objx;x++){   
            for(int y=1; y<=2*num_of_objy; y++){
                
                GLuint myLoc = glGetUniformLocation(gProgram, "randomnumber");
                glProgramUniform1i(gProgram, myLoc,objeler[x-1][y-1].color);
                
                float transy2=(float)-8+mr2+(y-2)*(2*mr2);
                if(objeler[x-1][y-1].buyuyor==true && objeler[x-1][y-1].scale<(m*1.5)){
                    float transx=(float)-10+mr+(x-1)*(2*mr);
                    float transy=(float)-8+mr2+(y-1)*(2*mr2);
                    glLoadIdentity();
                    glTranslatef(transx, transy, -10);
                    glRotatef(angle, 0, 1, 0);
                    glScalef(objeler[x-1][y-1].scale,objeler[x-1][y-1].scale,objeler[x-1][y-1].scale);
                    drawModel();
                    objeler[x-1][y-1].scale+=0.01;

                }
                else if(objeler[x-1][y-1].buyuyor==true && objeler[x-1][y-1].scale>=(m*1.5)){
                    
                    objeler[x-1][y-1].buyuyor=false;
                    objeler[x-1][y-1].boom=true;
                    objeler[x-1][y-1].scale =0;

                    
                    int y2=y+1;
                    
                    for(;y2<=2*num_of_objy; y2++){
                        objeler[x-1][y2-1].kayiyor=true;
                        
                    }
 
                    if(objeler[x-1][y-1].boom == true){
                        //objeler[x-1][y-1].scale = 0 ;
                    }
                    
                    drawModel();
                    
                    
                }
                else if(objeler[x-1][y-1].kayiyor==true && objeler[x-1][y-1].boom==false){
                    
                    float transx=(float)-10+mr+(x-1)*(2*mr);
                    float transy=(float)-8+mr2+(y-1)*(2*mr2);
                    
                    glLoadIdentity();
                    glTranslatef(transx, transy+objeler[x-1][y-1].sliding, -10);
                    glRotatef(angle, 0, 1, 0);
                    glScalef(m,m,m);
                    
                    objeler[x-1][y-1].sliding-=0.05;
                    
                    
                    drawModel();
                    
                    if(transy+objeler[x-1][y-1].sliding <=transy2) {
                        
                        objeler[x][y-1].color = objeler[x-1][y-1].color;
                        objeler[x-1][y-1].kayiyor=false;
                        objeler[x-1][y-2].boom=false;
                        //GLuint myLoc = glGetUniformLocation(gProgram, "randomnumber");
                        //glProgramUniform1i(gProgram, myLoc,objeler[x][y-1].color);
                        if(y>num_of_objy+1)
                        {   
                            objeler[x-1][y-1].boom=true;
                            
                        }
                        objeler[x-1][y-1].sliding=0;
                    }

                    
                }
                else{   
                    objeler[x-1][y-1].scale=m;
                        float transx=(float)-10+mr+(x-1)*(2*mr);
                        float transy=(float)-8+mr2+(y-1)*(2*mr2);
                        glLoadIdentity();
                        glTranslatef(transx, transy+ objeler[x-1][y-1].sliding, -10);
                        glRotatef(angle, 0, 1, 0);
                        glScalef(objeler[x-1][y-1].scale,objeler[x-1][y-1].scale,objeler[x-1][y-1].scale);
                        glColor3f(1,1,0);
                        drawModel();
                        
                } 
                
            }
            
        }
    //matchingFunction(a,b,m);
    drawModel();
    
	angle += 0.5;
}


void reshape(GLFWwindow* window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-10, 10, -10, 10, -20, 20);
    //gluPerspective(45, 1, 1, 100);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
}



void mainLoop(GLFWwindow* window, int a, int b)
{
    int toplam = 3*a*b;
    static float m=(float) (16*9)/(10*(a*b));
    for(int y = 0 ; y<3*a ; y++){
        vector<objeozellikleri> tempo3;
        for(int z=0; z<2*b;z++){
            int var_color = (rand()% 5);
            //std::cout<< var_color << endl;
            bool x=false;
            static float n=0;
            objeozellikleri tempo2(x,x,x,m,n,var_color);
            
            tempo3.push_back(tempo2);
        }
        objeler.push_back(tempo3);
    }



    while (!glfwWindowShouldClose(window))
    {
        
        if(clicked && mposx>=0 && mposx <=640 && mposy >= 0 && mposy<=600) {
            //buyumekte=true;
            objnumy= (mposy*b) / 540;
            objnumx=((mposx)*a)/640;
            objnumy=b-objnumy ;
            objnumy--;
            objeler[objnumx][objnumy].buyuyor=true;

        }
        display(a,b,m,toplam);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
 
}




int string_to_int(char *str){
    int i=0,x=0;
    for(;str[i];i++){
        if(str[i]>='0' && str[i]<='9'){
            x=x*10+(str[i]-48);
        }
        else 
            break;
    }
    return x;
}




//Ben Ekledim
static void cursorPositionCallback(GLFWwindow *window, double xPos, double yPos);
void mouseButtonCallback(GLFWwindow *window , int button, int action, int mods);


int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{   
    
    //cout << "Number of command line arguments (argc) entered: " << argc<<endl; 
    int a,b;
    a=string_to_int(argv[1]);
    b=string_to_int(argv[2]);
    string s=argv[3];/*
    if(s=="armadillo.obj"){

        for (int i = 0; i < argc; ++i) 
        cout <<"argv["<<i<<"] : "<<argv[1] << "\n"; 

    }*/
    
    //int n;
    //std::cin >> n;
    GLFWwindow* window;
    if (!glfwInit())
    {
        exit(-1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    //srand(time(0));
    int width = 640, height = 600;
    window = glfwCreateWindow(width, height, "Simple Example", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }
    glfwSetCursorPosCallback(window,cursorPositionCallback);//ben ekledim
    glfwSetMouseButtonCallback(window, mouseButtonCallback);

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    char rendererInfo[512] = {0};
    strcpy(rendererInfo, (const char*) glGetString(GL_RENDERER));
    strcat(rendererInfo, " - ");
    strcat(rendererInfo, (const char*) glGetString(GL_VERSION));
    glfwSetWindowTitle(window, rendererInfo);

    init(s);

    glfwSetKeyCallback(window, keyboard);
    glfwSetWindowSizeCallback(window, reshape);

    reshape(window, width, height); // need to call this once ourselves
    mainLoop(window,a,b); // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}

static void cursorPositionCallback(GLFWwindow *window, double xPos,double yPos) {
    mposx=xPos;
    mposy=yPos;

    if(xPos==200) std::cout<<mposx<<":"<<mposy<<std::endl;
    
}

void mouseButtonCallback(GLFWwindow *window , int button, int action, int mods){
    if(button==GLFW_MOUSE_BUTTON_LEFT && action==GLFW_PRESS){
        std::cout <<"Left button pressed"<<std::endl;

        clicked=true;
    }
    else{
        clicked=false;
    }
    if(clicked) std::cout<<mposx<<":"<<mposy<<std::endl;

}
