#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <list>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

bool visibltycontrol(double  denom, double num, double& d1, double& d2){
    double t = num;
    bool flag = true;
    if( !(denom >0) && !(denom < 0) && num > 0){
        flag = false;
    }
    else{
        if(denom > 0){
            if(t/denom > d2){
                flag = false;
            } else if(t/denom > d1){
                d1 = t/denom;
            }
        } 
        else if(denom < 0){
            if(t/denom < d1) flag = false;
            else if(t/denom < d2){
                d2 = t/denom;
            }
        } 
    }

    if(flag == false){
        return false;
    }
    return true;
}

int line_clipp(Vec4& ilkvector, Color& ilkrenk, Vec4& ikincivektor, Color& ikincirenk ) {
    
	int visibility = 0;
    double d1 = 0, d2 = 1;
    double xx = ikincivektor.x - ilkvector.x;
    double yy = ikincivektor.y - ilkvector.y;
    double zz = ikincivektor.z - ilkvector.z;
	Color cc = ikincirenk;
    cc= cc - ilkrenk;
    double minx , miny , minz ;
  
	bool forX_1 , forX_2 , forY_1 , forY_2 , forZ_1 , forZ_2;
	forX_1 = visibltycontrol(xx, -1-ilkvector.x, d1, d2);
	forX_2 = visibltycontrol(-xx, ilkvector.x-1, d1, d2);
	forY_1 = visibltycontrol(yy, -1-ilkvector.y, d1, d2);
	forY_2 = visibltycontrol(-yy, ilkvector.y-1, d1, d2);
	forZ_1 = visibltycontrol(zz, -1-ilkvector.z, d1, d2);
	forZ_2 = visibltycontrol(-zz, ilkvector.z-1, d1, d2);
    int calculating_number = forX_1 * forX_2 * forY_1 * forY_2 * forZ_1 * forZ_2;
    double calculating_numberFort =0;
    
    
    if(calculating_number){
        visibility = 1;
        if (d2 < 1){
        calculating_numberFort = d2;
        ikincivektor.x = ilkvector.x + (xx * calculating_numberFort);
        ikincivektor.y = ilkvector.y + (yy * calculating_numberFort);
        ikincivektor.z = ilkvector.z + (zz * calculating_numberFort);
        ikincirenk = ilkrenk + (cc * calculating_numberFort);
        
    }
        if (d1 > 0) {
        calculating_numberFort = d1;
        ilkvector.x = ilkvector.x + (xx * calculating_numberFort);
        ilkvector.y = ilkvector.y + (yy * calculating_numberFort);
        ilkvector.z = ilkvector.z + (zz * calculating_numberFort);
        ilkrenk = ilkrenk + (cc * calculating_numberFort);
    }
    
        
    }
            
   
    return visibility;
}



double minfinder(double array[3]){
    double min=array[0];
    int i;
    for (i = 1; i < 3; i++)
    {
        if (array[i] < min)
            min = array[i];
    }
    return min;
}

double maxfinder(double array[3]){
    double max=array[0];
    int i;
    for (i = 1; i < 3; i++)
    {
        if (array[i] >max)
            max = array[i];
    }
    return max;
}


void rasterize_triangle(vector<vector<Color>>& image, Vec4 vektor1,Vec4 vektor2,Vec4 vektor3,  Color  c0,  Color  c1, Color  c2, int height, int width){
    int min_x=0, min_y=0, max_x, max_y;
    double a,b,g;
    double array[3];
    array[0] = vektor1.x;
    array[1] = vektor2.x;
    array[2] = vektor3.x;

    Color color;
    Color color_helper;
    Color color_helper2;
    if(minfinder(array)>= 0){
        min_x = minfinder(array);
    }
    
    if(min_x > width -1){
        min_x = width -1 ;
    }

    array[0] = vektor1.y;
    array[1] = vektor2.y;
    array[2] = vektor3.y;
    if(minfinder(array) >= 0){
        min_y = minfinder(array);
    }

    if(min_y > height -1){
        min_y = height - 1;
    }

    array[0] = vektor1.x;
    array[1] = vektor2.x;
    array[2] = vektor3.x;



    max_x = maxfinder(array);
    if(max_x < 0) {
        max_x = 0;
    } 
    else if (max_x > width-1){
        max_x = width -1 ;

    }
    

    array[0] = vektor1.y;
    array[1] = vektor2.y;
    array[2] = vektor3.y;

    max_y = maxfinder(array);
    if(max_y < 0) {
        max_y = 0;
    } 
    else if (max_y > height-1){
        max_y = height -1 ;
    }


    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {
            a = ((x * (vektor2.y - vektor3.y)) + (y * (vektor3.x - vektor2.x)) + (vektor2.x * vektor3.y) - (vektor2.y * vektor3.x) )/ ((vektor1.x * (vektor2.y - vektor3.y)) + (vektor1.y * (vektor3.x - vektor2.x)) + (vektor2.x * vektor3.y) - (vektor2.y * vektor3.x));
            b = ((x * (vektor3.y - vektor1.y)) + (y * (vektor1.x - vektor3.x)) + (vektor3.x * vektor1.y) - (vektor3.y * vektor1.x)) / ((vektor2.x * (vektor3.y - vektor1.y)) + (vektor2.y * (vektor1.x - vektor3.x)) + (vektor3.x * vektor1.y) - (vektor3.y * vektor1.x));
            g = ((x * (vektor1.y - vektor2.y)) + (y * (vektor2.x - vektor1.x)) + (vektor1.x * vektor2.y) - (vektor1.y * vektor2.x) )/ ((vektor3.x * (vektor1.y - vektor2.y)) + (vektor3.y * (vektor2.x - vektor1.x)) + (vektor1.x * vektor2.y) - (vektor1.y * vektor2.x)) ;
            if (a >= 0 && b >= 0 && g >= 0){
                color = ((c0) * a);
                color_helper = ((c1) * b);
                color_helper2 = ((c2) * g);
                color_helper = color_helper + color_helper2;
                color = color + color_helper;
                image[x][y].r = (int)(color.r+0.5);
                image[x][y].g = (int)(color.g+0.5);
                image[x][y].b = (int)(color.b+0.5);
            }
        }
    }

}


void restarize_linehelper1( vector<vector<Color>> &img, Vec4 &vektorilk, Vec4 &vektorsecond, Color &colorilk, Color &coloriki , 
double &doublex , double &doubley , Color &renk1 , Color &color , int &d) {
        Vec4 temp;
        Color temp2;

        if (doublex<0) {
            
            temp = vektorilk;
            vektorilk = vektorsecond;
            vektorsecond = temp;
            temp2= colorilk;
            colorilk = coloriki;
            coloriki = temp2;
            doublex = vektorsecond.x - vektorilk.x;
            doubley = vektorsecond.y - vektorilk.y;
            
        }
        if (doubley<0) {
            
            int y = vektorilk.y;
            color = colorilk;
            d = ((-1)*doubley) - 0.5 * (doublex);
            renk1 = (coloriki - colorilk) / (doublex);
            for (int x = vektorilk.x; x <= vektorsecond.x; x++) {
                img[x][y].r = (int)(color.r+0.5);
                img[x][y].b = (int)(color.b+0.5);
                img[x][y].g = (int)(color.g+0.5);
                d = d+ ((-1)*doubley);
                if (d  > 0) { // choose NE
                    y-=1;
                    d = d - (doublex);
                }
                color = color + renk1;
            }
        }
        else {
            int y = vektorilk.y;
            color = colorilk;
            d = ((-1)*doubley) + (0.5 * (doublex));
            renk1 = (coloriki - colorilk) / (doublex);
            for (int x = vektorilk.x; x <= vektorsecond.x; x++) {
                img[x][y].r = (int)(color.r+0.5);
                img[x][y].b = (int)(color.b+0.5);
                img[x][y].g = (int)(color.g+0.5);
                d = d + ((-1)*doubley);
                if (d  < 0) {
                    y+=1;
                    d = d +  (doublex);
                } 
                    
                color = color + renk1;
            }
        }
}

void restarize_linehelper2( vector<vector<Color>> &img, Vec4 &vektorilk, Vec4 &vektorsecond, Color &colorilk, Color &coloriki , 
double &doublex , double &doubley , Color &renk1 , Color &color , int &d) {
        Vec4 temp;
        Color temp2;
        if (doubley<0) {
            temp = vektorilk;
            vektorilk = vektorsecond;
            vektorsecond = temp;
            // swap color
            temp2= colorilk;
            colorilk = coloriki;
            coloriki = temp2;
            doublex = vektorsecond.x - vektorilk.x;
            doubley = vektorsecond.y - vektorilk.y;
        }
        if (doublex<0) {
            int x = vektorilk.x;
            color = colorilk;
            d = (doublex) - 0.5 * ((-1)*doubley);
            renk1 = (coloriki - colorilk) / (doubley);

            for (int y = vektorilk.y; y <= vektorsecond.y; y++) {
                img[x][y].r = (int)(color.r+0.5);
                img[x][y].b = (int)(color.b+0.5);
                img[x][y].g = (int)(color.g+0.5);
                d = d + (doublex);
                if (d < 0) {
                    x --;
                    d = d -  ((-1)*doubley);
                }
                    
                color = color + renk1;
            }
        }
        else {
            int x = vektorilk.x;
            color = colorilk;
            d = (doublex) + ( -0.5 * (doubley));
            renk1 = (coloriki - colorilk) / (vektorsecond.y - vektorilk.y);

            for (int y = vektorilk.y; y <= vektorsecond.y; y++) {
                img[x][y].r = (int)(color.r+0.5);
                img[x][y].b = (int)(color.b+0.5);
                img[x][y].g = (int)(color.g+0.5);
                d = d + (doublex);
                if (d  > 0) {
                    x ++;
                    d =  d + (-1 * doubley);
                }
                    
                color = color + renk1;
            }
        }

}

Matrix4 Matrixer_translations(int i,Mesh *mesh,const vector<Translation *> translations, Matrix4 Matrixer){

    int id_s =(mesh->transformationIds)[i] - 1;
    Translation *translation = translations[id_s];
            double t_matrix[4][4] = {{1, 0, 0, translation->tx},
                                     {0, 1, 0, translation->ty},
                                     {0, 0, 1, translation->tz},
                                     {0, 0, 0, 1}};
            Matrix4 transMatrix(t_matrix);
            Matrixer = multiplyMatrixWithMatrix(transMatrix, Matrixer);
    return Matrixer;

}
Matrix4 Matrixer_scalings(int i,Mesh *mesh,const vector<Scaling *> scalings, Matrix4 Matrixer){

        int id_s = (mesh->transformationIds[i]) - 1;
        Scaling *scaling = scalings[id_s];
        double s_matrix[4][4] = {{scaling->sx, 0,0,0},
                                     {0,scaling->sy, 0,0},
                                     {0,0,scaling->sz, 0},
                                     {0,0,0,1}};
            Matrix4 scalingMatrix(s_matrix);
            Matrixer = multiplyMatrixWithMatrix(scalingMatrix, Matrixer);
        return Matrixer;
}


Matrix4 Matrixer_rotations_helper(Vec3 &u ,Vec3 &w , Vec3 &v, Rotation * rotation,Matrix4 Matrixer){

    double _Matrix1[4][4] = {{u.x, u.y, u.z, 0},
                              {v.x, v.y, v.z, 0},
                              {w.x, w.y, w.z, 0},
                              {0 ,0,0,1}};

    Matrix4 _Matrix2(_Matrix1);

    double _Matrix3[4][4] = {{u.x, v.x, w.x, 0},
                                  {u.y, v.y, w.y, 0},
                                  {u.z, v.z, w.z, 0},
                                  {0,   0,   0,   1}};
        

    Matrix4 _Matrix4(_Matrix3);

    double theta = (rotation->angle * M_PI) / (180.0);
        
    double _Matrix5[4][4] = {{1, 0,   0,  0},
                        {0, cos(theta), (-1) * sin(theta), 0},
                        {0, sin(theta), cos(theta),   0},
                        {0, 0, 0,  1}};

    Matrix4 _Matrix6(_Matrix5);
        

        Matrix4 m1 = multiplyMatrixWithMatrix(_Matrix6, _Matrix2);
        Matrix4 m2 = multiplyMatrixWithMatrix(_Matrix4, m1);
        Matrix4 helper_Matrix2 = multiplyMatrixWithMatrix(m2, Matrixer);

        return helper_Matrix2;

}
double CalculateAbsoluteValue(double value){

    if(value<0){
        value = value *(-1);
    }
    return value;
}

Matrix4 Matrixer_rotations(int i,Mesh *mesh,const vector<Rotation *> rotations, Matrix4 Matrixer){
     
        int id_s =mesh->transformationIds[i] - 1;
        Rotation *rotation = rotations[id_s];
        Vec3 u, v, new_Vec;
        u.x = rotation->ux;
        u.y = rotation->uy;
        u.z = rotation->uz;
        u.colorId = -1;
        double absolute_ux = CalculateAbsoluteValue(u.x);
        double absolute_uy = CalculateAbsoluteValue(u.y);
        double absolute_uz = CalculateAbsoluteValue(u.z);
        double new_array[3];
        new_array[0] = absolute_ux;
        new_array[1] = absolute_uy;
        new_array[2] = absolute_uz;
        double min_val = minfinder(new_array);
        if (min_val == absolute_ux) {
            v.x = 0.0;
            v.y = (-1.0) * u.z;
            v.z = u.y;
            v.colorId = -1;
        } 
        else if (min_val == absolute_uy) {
            v.x = (-1.0) * u.z;
            v.y = 0.0;
            v.z = u.x;
            v.colorId = -1;
        } 
        else if (min_val == absolute_uz) {
            v.x = (-1.0) * u.y;   
            v.y = u.x;
            v.z = 0.0;
            v.colorId = -1;
        } 
        
        new_Vec = crossProductVec3(u, v);
        new_Vec = normalizeVec3(new_Vec);
        v = normalizeVec3(v);

        Matrixer = Matrixer_rotations_helper(u , new_Vec ,v,rotations[id_s] , Matrixer);

        return Matrixer;
}


Matrix4 Matrixer(Mesh *mesh, const vector<Scaling *> scalings,const vector<Rotation *> rotations,const vector<Translation *> translations) {

    Matrix4 Matrixer = getIdentityMatrix();
    for (int i = 0; i < mesh->numberOfTransformations; i++) {
        if (((mesh->transformationTypes)[i]) == 't') {
           Matrixer = Matrixer_translations(i,mesh, translations,  Matrixer);
        } 
        else if (((mesh->transformationTypes)[i]) == 's') {
            Matrixer = Matrixer_scalings(i,mesh, scalings,  Matrixer); 
        } 
        else if (((mesh->transformationTypes)[i]) == 'r') {
            Matrixer = Matrixer_rotations(i,mesh, rotations,  Matrixer);
        }
    }

    return Matrixer;
}



Matrix4 pipelineHelperCam(Camera *camera){
    
    double TransMatrix[4][4] =
                {{1 , 0, 0, -1*(camera -> pos.x)},
                 {0,1,0, -1*(camera -> pos.y)},
                 {0,0,1,-1*(camera -> pos.z)},
                 {0,0,0,1}};
    double RotationMatrix[4][4] =
                {{camera -> u.x, camera -> u.y, camera -> u.z,0},
                 {camera -> v.x, camera -> v.y, camera -> v.z,0},
                 {camera -> w.x, camera -> w.y, camera -> w.z ,0},
                 {0,0,0,1} };

    return (multiplyMatrixWithMatrix(RotationMatrix, TransMatrix));
}


Matrix4 pipelineHelperProject(Camera *camera){
    if(camera -> projectionType){
        double PerspectiveM[4][4] =
                {{((2*camera -> near) / (camera -> right - camera -> left)), 0, ((camera -> right + camera -> left) / (camera -> right - camera -> left)),0},
                 {0, ((2 * camera -> near) / (camera -> top - camera -> bottom)), ((camera -> top + camera -> bottom)/(camera -> top - camera -> bottom)), 0},
                 {0,0,-1*((camera -> far + camera -> near) / (camera -> far - camera -> near)), -1*((2*camera->far *  camera->near) / (camera->far - camera->near))},
                 {0,0,-1, 0}};
        return (Matrix4(PerspectiveM));
    } else {
    double  OrthogonalM[4][4] =
            {{2/(camera -> right - camera->left), 0, 0, -1*((camera -> right + camera -> left) / (camera -> right - camera -> left))},
             {0, 2/(camera -> top - camera -> bottom), 0, -1*((camera -> top + camera -> bottom) / (camera -> top - camera -> bottom))},
             {0, 0, -1*(2/(camera -> far - camera -> near)), -1*((camera -> far + camera -> near) / (camera -> far - camera -> near))},
             {0, 0, 0, 1}};
        return (Matrix4(OrthogonalM));
    }
}


Matrix4 pipelineHelperView(Camera *camera){
    double ViewportM[4][4]=
            {{(camera -> horRes / (double)2),0,0,((camera -> horRes -1)/(double)2) },
             {0,(camera -> verRes/(double)2), 0, ((camera -> verRes -1)/(double)2) },
             {0,0,0.5,0.5},
             {0, 0, 0, 1}};

    return (Matrix4(ViewportM));
}




void Scene::forwardRenderingPipeline(Camera *camera)
{
    Matrix4 cameraninmatrixi;
    cameraninmatrixi = pipelineHelperCam(camera);

    Matrix4 projectionmatrixi;
    projectionmatrixi = pipelineHelperProject(camera);

    Matrix4 viewmatrixi;
    viewmatrixi = pipelineHelperView(camera) ;


    for (int i=0;i<this->meshes.size();i++) {
        Mesh *mesh=meshes[i];
        Matrix4 sonmatrixi = multiplyMatrixWithMatrix(projectionmatrixi, multiplyMatrixWithMatrix(cameraninmatrixi, Matrixer(mesh,  this->scalings, this->rotations,this->translations))); // Mproj * Mcam * Mmodel
        for (int j = 0; j < mesh->numberOfTriangles; j++) {
            const Vec3 *vect0 = this->vertices[((mesh->triangles[j]).getFirstVertexId()) - 1];
            const Vec3 *vect1 = this->vertices[((mesh->triangles[j]).getSecondVertexId()) - 1];
            const Vec3 *vect2 = this->vertices[((mesh->triangles[j]).getThirdVertexId()) - 1];

            const Color *col0 = this->colorsOfVertices[vect0->colorId - 1];
            const Color *col1 = this->colorsOfVertices[vect1->colorId - 1];
            const Color *col2 = this->colorsOfVertices[vect2->colorId - 1];
            

            Vec4 firstpromatrixi = multiplyMatrixWithVec4(sonmatrixi, Vec4(vect0->x, vect0->y, vect0->z, 1, vect0->colorId));
            Vec4 ikincipromatrixi = multiplyMatrixWithVec4(sonmatrixi, Vec4(vect1->x, vect1->y, vect1->z, 1, vect1->colorId));
            Vec4 ucuncupromatrixi = multiplyMatrixWithVec4(sonmatrixi, Vec4(vect2->x, vect2->y, vect2->z, 1, vect2->colorId));

                bool flag =true;
                Vec3 vv_0, vv_1, vv_2;
                Vec4* mylist = new Vec4[3];
                mylist[0] = firstpromatrixi;
                mylist[1] = ikincipromatrixi;
                mylist[2] = ucuncupromatrixi;

                Vec3* mylist2 = new Vec3[3];
                mylist2[0] = vv_0;
                mylist2[1] = vv_1;
                mylist2[2] = vv_2;

                for(int i = 0 ; i<3 ; i++){
                    mylist2[i].x = mylist[i].x;
                    mylist2[i].y = mylist[i].y;
                    mylist2[i].z = mylist[i].z;
                    mylist2[i].colorId = mylist[i].colorId;
                }

                    Vec3 normalizedvec = normalizeVec3(crossProductVec3(subtractVec3(mylist2[1], mylist2[0]), subtractVec3(mylist2[2], mylist2[0])));
                    double mrx = dotProductVec3(normalizedvec,mylist2[0]);
                    if(mrx<0){
                        flag=true;
                    }
                    else{
                        flag = false;
                    }

            if (this->cullingEnabled && flag) {
                continue;
            }
            firstpromatrixi.dividebyt();
            ikincipromatrixi.dividebyt();
            ucuncupromatrixi.dividebyt();
            
            if (!mesh->type) {

                Color colorfirstsecond_0= *col0 , colorfirstsecond_1 = *col1 , colorthirdfirst_0= *col2 , colorthirdfirst_1 = *col0 , colorsecondthird_0= *col1 , colorsecondthird_1 = *col2;
                Vec4 vectorfirstsecon_0  ,vectorfirstsecond_1  , vectorsecondthird_0  ,vectorsecondthird_1 ,  vectorthirdfirst_0  ,vectorthirdfirst_1  ;

                vectorfirstsecon_0 = firstpromatrixi;
                vectorthirdfirst_1 = vectorfirstsecon_0;

                vectorfirstsecond_1 =  ikincipromatrixi;
                vectorsecondthird_0 = vectorfirstsecond_1;

                vectorsecondthird_1 = ucuncupromatrixi;
                vectorthirdfirst_0 = vectorsecondthird_1;
                int tempo1 =line_clipp(vectorfirstsecon_0, colorfirstsecond_0, vectorfirstsecond_1, colorfirstsecond_1);
                int tempo2 = line_clipp(vectorsecondthird_0, colorsecondthird_0, vectorsecondthird_1, colorsecondthird_1);
                int tempo3 = line_clipp(vectorthirdfirst_0, colorthirdfirst_0, vectorthirdfirst_1, colorthirdfirst_1);
                bool mr1=false ,mr2=false ,mr3=false ;
                if(tempo1==1) mr1=true;
                if(tempo2==1) mr2=true;
                if(tempo3==1) mr3=true;

                

                if (mr1) {
                    vectorfirstsecon_0 = multiplyMatrixWithVec4(viewmatrixi, vectorfirstsecon_0);
                    vectorfirstsecond_1 = multiplyMatrixWithVec4(viewmatrixi, vectorfirstsecond_1);
                    vectorsecondthird_0 = multiplyMatrixWithVec4(viewmatrixi, vectorsecondthird_0);
                    vectorsecondthird_1 = multiplyMatrixWithVec4(viewmatrixi, vectorsecondthird_1);
                    vectorthirdfirst_0 = multiplyMatrixWithVec4(viewmatrixi, vectorthirdfirst_0);
                    vectorthirdfirst_1 = multiplyMatrixWithVec4(viewmatrixi, vectorthirdfirst_1);

                    double doublex = vectorfirstsecond_1.x - vectorfirstsecon_0.x;
                    double doubley = vectorfirstsecond_1.y - vectorfirstsecon_0.y;
                    double m = doubley/ doublex;
                    
                    double x1=0, y1, z1, t1;
                    
                    Color temp2;
                    int d;
                    Color renk1, color;

                    if ((m>=-1 && m<=1)) {
                        restarize_linehelper1( this->image, vectorfirstsecon_0, vectorfirstsecond_1, colorfirstsecond_0, colorfirstsecond_1 , 
                doublex , doubley , renk1 , color , d);
                    }
                    else {
                        restarize_linehelper2( this->image, vectorfirstsecon_0, vectorfirstsecond_1, colorfirstsecond_0, colorfirstsecond_1 , 
                doublex , doubley , renk1 , color , d);
                    }
                }
                else {
                    vectorfirstsecon_0 = multiplyMatrixWithVec4(viewmatrixi, vectorfirstsecon_0);
                vectorfirstsecond_1 = multiplyMatrixWithVec4(viewmatrixi, vectorfirstsecond_1);
                vectorsecondthird_0 = multiplyMatrixWithVec4(viewmatrixi, vectorsecondthird_0);
                vectorsecondthird_1 = multiplyMatrixWithVec4(viewmatrixi, vectorsecondthird_1);
                vectorthirdfirst_0 = multiplyMatrixWithVec4(viewmatrixi, vectorthirdfirst_0);
                vectorthirdfirst_1 = multiplyMatrixWithVec4(viewmatrixi, vectorthirdfirst_1);
                }
                if (mr2) {

                    double doublex = vectorsecondthird_1.x - vectorsecondthird_0.x;
                    double doubley = vectorsecondthird_1.y - vectorsecondthird_0.y;
                    double mm = doubley/ doublex;
                    
                    double x1=0, y1, z1, t1;
                    
                    Color temp2;
                    int d;
                    Color renk1, color;

                    if ((mm>=-1 && mm<=1)) {
                        restarize_linehelper1( this->image, vectorsecondthird_0, vectorsecondthird_1, colorsecondthird_0, colorsecondthird_1 , 
                doublex , doubley , renk1 , color , d);
                    }
                    else {
                        restarize_linehelper2( this->image, vectorsecondthird_0, vectorsecondthird_1, colorsecondthird_0, colorsecondthird_1 , 
                doublex , doubley , renk1 , color , d);
                    }
                }

                if (mr3) {

                   double doublex = vectorthirdfirst_1.x - vectorthirdfirst_0.x;
                    double doubley = vectorthirdfirst_1.y - vectorthirdfirst_0.y;
                    double mm = doubley/ doublex;
                    
                    double x1=0, y1, z1, t1;
                    
                    Color temp2;
                    int d;
                    Color renk1, color;

                    if ((mm>=-1 && mm<=1)) {
                        restarize_linehelper1( this->image, vectorthirdfirst_0, vectorthirdfirst_1, colorthirdfirst_0, colorthirdfirst_1 , 
                doublex , doubley , renk1 , color , d);
                    }
                    else {
                        restarize_linehelper2( this->image, vectorthirdfirst_0, vectorthirdfirst_1, colorthirdfirst_0, colorthirdfirst_1 , 
                doublex , doubley , renk1 , color , d);
                    }


                }

            } 
            else {

                rasterize_triangle( this->image, multiplyMatrixWithVec4(viewmatrixi, firstpromatrixi), multiplyMatrixWithVec4(viewmatrixi, ikincipromatrixi), multiplyMatrixWithVec4(viewmatrixi, ucuncupromatrixi),  *col0, *col1, *col2, camera->verRes, camera->horRes);

            }
        }
    }

}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL) {
		str = pElement->GetText();
		
		if (strcmp(str, "enabled") == 0) {
			cullingEnabled = true;
		}
		else {
			cullingEnabled = false;
		}
	}

	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		// read projection type
		str = pCamera->Attribute("type");

		if (strcmp(str, "orthographic") == 0) {
			cam->projectionType = 0;
		}
		else {
			cam->projectionType = 1;
		}

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read meshes
	pElement = pRoot->FirstChildElement("Meshes");

	XMLElement *pMesh = pElement->FirstChildElement("Mesh");
	XMLElement *meshElement;
	while (pMesh != NULL)
	{
		Mesh *mesh = new Mesh();

		pMesh->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = pMesh->Attribute("type");

		if (strcmp(str, "wireframe") == 0) {
			mesh->type = 0;
		}
		else {
			mesh->type = 1;
		}

		// read mesh transformations
		XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *clone_str;
		int v1, v2, v3;
		XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
		clone_str = strdup(str);

		row = strtok(clone_str, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);
			
			if (result != EOF) {
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		meshes.push_back(mesh);

		pMesh = pMesh->NextSiblingElement("Mesh");
	}
}

/*

	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}


