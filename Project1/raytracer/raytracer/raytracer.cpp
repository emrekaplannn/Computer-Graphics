#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "helper_functions.h"

typedef unsigned char RGB[3];

class Camera {
public:
    int _width, _height; // resolution of camera // nx = width
    float _distance; // distance = near_distance
    Vector_3d _position; // position of camera
    Vector_3d gaze;
    parser::Vec4f nearPlane; // x = left , y = right , z = top , w = bottom
    Vector_3d up_3d, up_u, up_w; //up_3d = v , up_u = u , up_w = w
    std::string image_name;

    Ray generateRay(int i, int j, Vector_3d bgcolor);
};

class Ray {
public:
    
    Vector_3d _origin, _direction, color,bg;
    
    Vector_3d claculatediffuse(Vector_3d& diffuseprovider,  Vector_3d& M, const Vector_3d& R, Vector_3d &irradiance);
    Vector_3d calculatespecular(float phong_exp,Vector_3d &spec_coef , Vector_3d &yarim, const Vector_3d &MR,Vector_3d &irradiance);
    
    void computecolour(Mesh *meshes,int num_meshes,Triangle *triangles,int num_triangles,Sphere *spheres,int num_spheres,Light *lights,int num_lights,int max_recursion,float epsilon);

    void computeColorForSphere(Mesh *meshes, int num_meshes, Triangle *triangles, int num_triangles, Sphere *s, int num_spheres,
    Light *light, int num_light, float epsilon , float minT, Vector_3d &L , Vector_3d N, Vector_3d P , int minI);
    
    void computeColorForTriangle(Mesh *meshes, int num_meshes, Triangle *triangles, int num_triangles, Sphere *s, int num_spheres,
    Light *light, int num_light, float epsilon ,  float minT,  Vector_3d &L , Vector_3d N, Vector_3d P , int minI);
    
    void computeColorForMesh(Mesh *meshes, int num_meshes, Triangle *triangles, int num_triangles, Sphere *s, int num_spheres,
    Light *light, int num_light,float epsilon ,  float minT,  Vector_3d &L , Vector_3d N, Vector_3d P , int minI);

};

Ray
Camera::generateRay(int i, int j,Vector_3d bgcolor) {
    Ray result;
    float su, sv;
    Vector_3d m, a;
    su = (i + 0.5) * (this->nearPlane.y - this->nearPlane.x)/ _width;
    sv = (j + 0.5) * (this->nearPlane.z - this->nearPlane.w)/ _height;
    m = this->_position + gaze * _distance ;
    a = m + ((this->up_u) * (this->nearPlane.x) + (this->up_3d) *(this->nearPlane.z)) + ((this->up_u) * (su) + (this->up_3d)* (sv * -1));
    result._origin = this->_position;
    result._direction = a + _position * (-1);
    result.bg=bgcolor;
    return result;

}
class Mesh {
public:
    Vector_3d ambient_light,normal , ambient_ref,diffuse_ref,spec_ref,mirror_ref;
    float mr_exp;
    std::vector<Triangle> triangles;
    float intersect(Ray ray);
};

class Light {
public:
    Vector_3d position,_instensity;
    Vector_3d irradiance(const Vector_3d &_point){
        Vector_3d result;
        float _distance = sqrt((position.x - _point.x) * (position.x - _point.x) + (position.y - _point.y) * (position.y - _point.y) + (position.z - _point.z) * (position.z - _point.z));
        _distance = _distance * _distance;
        if (_distance != 0) {
        result.x = (this->_instensity.x) / _distance;
        result.y = (this->_instensity.y) / _distance;
        result.z = (this->_instensity.z) / _distance;
        }
        return result;
    }
};

class Triangle {
public:
    Vector_3d ambient_light ,ambient_ref,diffuse_ref,spec_ref,mirror_ref , first,second,third,normal;
    int tri;
    float mr_exp;
    float intersect(Ray ray);
};


float Triangle::intersect(Ray ray) {
    Vector_3d mr1((first.x - second.x), (first.y - second.y), (first.z - second.z)),
            mr2((first.x - third.x), (first.y - third.y), (first.z - third.z));
    float detfirst = mr1.x * ((mr2.y * ray._direction.z) - (ray._direction.y * mr2.z)) - mr2.x * ((mr1.y * ray._direction.z) - (ray._direction.y * mr1.z)) + ray._direction.x * ((mr1.y * mr2.z) - (mr2.y * mr1.z));
    Vector_3d mr11((first.x - ray._origin.x), (first.y - ray._origin.y), (first.z - ray._origin.z)),
            mr12((first.x - third.x), (first.y - third.y), (first.z - third.z));
    float detsecond = mr11.x * ((mr12.y * ray._direction.z) - (ray._direction.y * mr12.z)) - mr12.x * ((mr11.y * ray._direction.z) - (ray._direction.y * mr11.z)) + ray._direction.x * ((mr11.y * mr12.z) - (mr12.y * mr11.z));
    float detmr = mr1.x * ((mr11.y * ray._direction.z) - (ray._direction.y * mr11.z)) - mr11.x * ((mr1.y * ray._direction.z) - (ray._direction.y * mr1.z)) + ray._direction.x * ((mr1.y * mr11.z) - (mr11.y * mr1.z));
    float detthird = mr1.x * ((mr12.y * mr11.z) - (mr11.y * mr12.z)) - mr12.x * ((mr1.y * mr11.z) - (mr11.y * mr1.z)) + mr11.x * ((mr1.y * mr12.z) - (mr12.y * mr1.z));
    
    if (detfirst == 0) return -1;
    detsecond /= detfirst;
    detmr /= detfirst;
    detthird /= detfirst;

    if (detsecond >= -0.001 && detsecond + detmr <= 1 && detmr >= -0.001 ) return detthird;    
    else return -1;

}

Vector_3d Ray::claculatediffuse(Vector_3d &diffuseprovider,  Vector_3d &M, const Vector_3d &R,
                                Vector_3d &irradiance)   {             
    Vector_3d result = irradiance * 0;
    if (M * R > 0) result = irradiance * (M*R);
    result.x *= diffuseprovider.x;
    result.y *= diffuseprovider.y;
    result.z *= diffuseprovider.z;
    return result;
}

Vector_3d Ray::calculatespecular(float phong_exp, Vector_3d &spec_coef, Vector_3d &yarim,
                                 const Vector_3d &MR, Vector_3d &irradiance) {
    float phong = yarim * MR;
    if (phong < 0) phong = 0;
    phong = pow(phong, phong_exp);
    Vector_3d result = spec_coef * phong;
    result.x = result.x * irradiance.x;
    result.y = result.y * irradiance.y;
    result.z = result.z * irradiance.z;

    return result;
}


class Sphere{
public:
    Vector_3d ambient_light , center, color , ambient_ref,diffuse_ref,spec_ref,mirror_ref;
    float mr_exp , radius , mrconst;
    float intersect(Ray ray);
};

float Sphere::intersect(Ray ray) { //tolga hoca ile birebir aynı
    float A, B, C;
    float delta;
    float t, t1, t2;
    C = (ray._origin.x - center.x) * (ray._origin.x - center.x) + (ray._origin.y - center.y) * (ray._origin.y - center.y) +
        (ray._origin.z - center.z) * (ray._origin.z - center.z) - radius * radius;

    B = 2 * (ray._direction.x) * (ray._origin.x - center.x) + (2 * ray._direction.y) * (ray._origin.y - center.y) +
        (2 * ray._direction.z) * (ray._origin.z - center.z);
    A = ray._direction.x * ray._direction.x + ray._direction.y * ray._direction.y + ray._direction.z * ray._direction.z;
    delta = B * B - 4 * A * C;

    if (delta < 0) return -1;

    else if (delta == 0) {
        t = -B / (2 * A);
    } else {
        delta = sqrt(delta);
        A = 2 * A;
        t1 = (-B + delta) / A;
        t2 = (-B - delta) / A;
        if (t1 < t2) t = t1; else t = t2;
    }
    return t;
}
// for mesh
float Mesh::intersect(Ray ray) {
    float min = 90000, t ,x;
    int booll = 0;
    for (int i = 0; i < this->triangles.size(); i++) {
        t = triangles[i].intersect(ray);

        if (t < min && t >= 0) {
            min = t;
            booll = 1;
            x=1;
            this->normal = triangles[i].normal;
        }
    }
    if (booll)   return min;
    else    return -1;
}
void Ray::computeColorForSphere(Mesh *meshes, int num_meshes, Triangle *triangles, int num_triangles, Sphere *s, int num_spheres,
Light *light, int num_light, float epsilon , float value, Vector_3d &L , Vector_3d N, Vector_3d P , int firstTouch){

        color.x= s[firstTouch].ambient_ref.x * s[firstTouch].ambient_light.x;
        color.y= s[firstTouch].ambient_ref.y * s[firstTouch].ambient_light.y;
        color.z= s[firstTouch].ambient_ref.z * s[firstTouch].ambient_light.z;
        P = this->_origin + ((this->_direction) * value);
        N = (P + (s[firstTouch].center) * -1).normalize();

        for (int i = 0; i < num_light; i++) {
            
            float shadow_value = -1, light_value;
            Vector_3d light_3D_vector;
            L = light[i].position + P * (-1);
            light_value = sqrt(L*L);
            L = L.normalize();
            bool selector = false;

            Ray forshadow;
            forshadow._direction = L;
            forshadow._origin = P + N * (epsilon);
            light_3D_vector = light[i].position + P * (-1);

            for (int j = 0; j < num_spheres; j++) {

                if (j != firstTouch) {
                    shadow_value = s[j].intersect(forshadow);

                }
                if (shadow_value != -1) {
                    selector = true;
                    goto spherehere;
                }
            }
            spherehere:
            if (selector == false) {

                for (int j = 0; j < num_meshes; j++) {
                    shadow_value = meshes[j].intersect(forshadow);
                    if (shadow_value != -1) {
                        selector = true;
                        goto spherehere2;
                    }
                }
            }
            spherehere2:
            if (selector == false) {
                for (int j = 0; j < num_triangles; j++) {
                    shadow_value = triangles[j].intersect(forshadow);
                    if (shadow_value != -1 ) {
                        selector = true;
                        goto spherehere3;
                    }
                }
            }
            spherehere3:

            if (shadow_value <= light_value && shadow_value != -1 && shadow_value >= 0 && selector == true) {
                //do nothing
            } 
            else{
                Vector_3d irradiance = light[i].irradiance(P);
                Vector_3d rayd = L + ((this->_direction).normalize() * (-1));
                Vector_3d norm = rayd.normalize();
                Vector_3d diffuse = claculatediffuse(s[firstTouch].diffuse_ref, L, N, irradiance);
                Vector_3d specular = calculatespecular(s[firstTouch].mr_exp, s[firstTouch].spec_ref, norm, N, irradiance);
                this->color = this->color + diffuse;
                this->color = this->color + specular;
            }
        }


}
void Ray::computeColorForTriangle(Mesh *meshes, int num_meshes, Triangle *triangles, int num_triangles, Sphere *s, int num_spheres,
    Light *light, int num_light, float epsilon , float value, Vector_3d &L , Vector_3d N, Vector_3d P , int firstTouch){
        
        color.x= triangles[firstTouch].ambient_ref.x * triangles[firstTouch].ambient_light.x;
        color.y= triangles[firstTouch].ambient_ref.y * triangles[firstTouch].ambient_light.y;
        color.z= triangles[firstTouch].ambient_ref.z * triangles[firstTouch].ambient_light.z;
        P = this->_origin + (this->_direction) * value;
        N = N.normalize();
        for (int i = 0; i < num_light; i++) {
            float shadow_value = -1, light_value;
            L = light[i].position + P * (-1);
            light_value = sqrt(L*L);
            L = L.normalize();
            bool selector = false;
            Ray forshadow;
            forshadow._direction = L;
            forshadow._origin = P + N * (epsilon);
            for (int j = 0; j < num_spheres; j++) {
                shadow_value = s[j].intersect(forshadow);
                if (shadow_value != -1) {
                    selector = true;
                    goto trianglehere;
                }
            }
            trianglehere:
            if (selector==false) {
                for (int j = 0; j < num_meshes; j++) {

                    shadow_value = meshes[j].intersect(forshadow);

                    if (shadow_value != -1) {
                        selector = true;
                        goto trianglehere2;
                    }
                }
            }
            trianglehere2:
            if (selector==false) {
                for (int j = 0; j < num_triangles; j++) {
                    if (j != firstTouch)
                        shadow_value = triangles[j].intersect(forshadow);
                    if (shadow_value != -1) {
                        selector = true;
                        goto trianglehere3;
                    }
                }
            }
            trianglehere3:

            if (shadow_value <= light_value && shadow_value != -1 && shadow_value >= 0 && selector == true) {
                //do nothing
            } 
            else{
                Vector_3d irradiance = light[i].irradiance(P);
                Vector_3d rayd = L + ((this->_direction).normalize() * (-1));
                rayd = rayd.normalize();
                Vector_3d diffuse = claculatediffuse(triangles[firstTouch].diffuse_ref, L, N, irradiance);
                Vector_3d specular = calculatespecular(triangles[firstTouch].mr_exp, triangles[firstTouch].spec_ref, rayd, N, irradiance);
                this->color = this->color + diffuse;
                this->color = this->color + specular;
            }
        }
    }
    
void Ray::computeColorForMesh(Mesh *meshes, int num_meshes, Triangle *triangles, int num_triangles, Sphere *s, int num_spheres,
    Light *light, int num_light, float epsilon , float value, Vector_3d &L , Vector_3d N, Vector_3d P , int firstTouch){
        
        color.x= meshes[firstTouch].ambient_ref.x * meshes[firstTouch].ambient_light.x;
        color.y= meshes[firstTouch].ambient_ref.y * meshes[firstTouch].ambient_light.y;
        color.z= meshes[firstTouch].ambient_ref.z * meshes[firstTouch].ambient_light.z;
        N = meshes[firstTouch].normal;
        P = this->_origin + ((this->_direction)) * value;
        N = N.normalize();

        for (int i = 0; i < num_light; i++) {
            float shadow_value = -1, light_value;

            L = light[i].position + P * (-1);
            light_value = sqrt(L*L);
            L = L.normalize();


            bool selector = false;
            Ray forshadow;
            forshadow._direction = L;
            forshadow._origin = P + N * (epsilon);
            for (int j = 0; j < num_spheres; j++) {
                shadow_value = s[j].intersect(forshadow);
                if (shadow_value != -1 && shadow_value > 0) {
                    selector = true;
                    goto meshhere;
                }
            }
            meshhere:
            if (selector == false) {
                for (int j = 0; j < num_meshes; j++) {

                    shadow_value = meshes[j].intersect(forshadow);

                    if (shadow_value != -1) {
                        selector = true;
                        goto meshhere2;
                    }
                }
            }
            meshhere2:
            if (selector == false) {
                for (int j = 0; j < num_triangles; j++) {
                    shadow_value = triangles[j].intersect(forshadow);
                    if (shadow_value != -1 && shadow_value > 0)        
                    {
                        selector = true;
                        goto meshhere3;
                    }
                }
            }
            meshhere3:
                
            if (shadow_value < light_value && shadow_value != -1 && shadow_value >= 0 && selector == true) {
                //do nothing
            }
            else{
                Vector_3d irradiance = light[i].irradiance(P);
                Vector_3d rayd = L + ((this->_direction).normalize() * (-1));
                rayd = rayd.normalize();
                Vector_3d diffuse = claculatediffuse(meshes[firstTouch].diffuse_ref, L, N, irradiance);
                Vector_3d specular = calculatespecular(meshes[firstTouch].mr_exp, meshes[firstTouch].spec_ref, rayd, N, irradiance);
                this->color = this->color + diffuse;
                this->color = this->color + specular;
            }    
        }
    }

//  pixels
void Ray::computecolour(Mesh *meshes, int num_meshes, Triangle *triangles, int num_triangles, Sphere *s, int num_spheres,
                       Light *light, int num_light, int max_recursion,float epsilon) {
    
    float minT = 90000;
    int material_id_selector;
    bool alert = false;
    float distance_value;
    Vector_3d L, N, P;
    int minI;
    minI = -1;
    int minI_formirror = -1;
    int minT_formirror = 90000;
    int material_id_selector_formirror = 0;
    this->color = this->bg;
    for (int i = 0; i < num_spheres; i++) {
        distance_value = s[i].intersect(*this);
        if (distance_value < minT && distance_value > 0) {
            material_id_selector = 1;
            minT = distance_value;
            minI = i;
            alert=true;
        }
    }
    if (alert) {
        if(minI != -1){
            computeColorForSphere(meshes, num_meshes, triangles, num_triangles, s, num_spheres,
                        light, num_light,  epsilon ,  minT, L , N, P , minI);

            P = this->_origin + (this->_direction) * minT;

            N = P + (s[minI].center) * -1;
            N = N.normalize();
        }
    }
    alert = false;
    for (int i = 0; i < num_triangles; i++) {
        N = triangles[i].normal;
        distance_value = (triangles[i]).intersect(*this);
        if (distance_value < minT && distance_value > 0) {
            material_id_selector = 2;
            minT = distance_value;
            minI = i;
            alert=true;
        }
    }
    if (alert) {
        if(minI != -1){
            computeColorForTriangle(meshes, num_meshes, triangles, num_triangles, s, num_spheres,
                        light, num_light,  epsilon , minT,  L , N, P , minI);
                        P = this->_origin + (this->_direction) * minT;
            N = N.normalize();
        }
    }
    alert =false;
    for (int i = 0; i < num_meshes; i++) {
        distance_value = (meshes[i]).intersect(*this);
        if (distance_value < minT && distance_value > 0) {
            material_id_selector = 3;
            minT = distance_value;
            minI = i;
            alert=true;
        }
    }
    if (alert){
        if(minI != -1){
            computeColorForMesh(meshes, num_meshes, triangles, num_triangles, s, num_spheres,
                        light, num_light,  epsilon ,  minT,  L , N, P , minI);
                        N = meshes[minI].normal;
            P = this->_origin + ((this->_direction)) * minT;
            N = N.normalize();
        }
    }
    
    
    bool mirror_existed = 0;
    if (material_id_selector == 1 && minI != -1){
        if(s[minI].mirror_ref.x > 0 || s[minI].mirror_ref.y > 0 || s[minI].mirror_ref.z > 0)  mirror_existed = 1;

    }
    if (material_id_selector == 2 && minI != -1){
        if (triangles[minI].mirror_ref.x > 0 || triangles[minI].mirror_ref.y > 0 || triangles[minI].mirror_ref.z > 0) mirror_existed = 1;
        
    }
    if (material_id_selector == 3 && minI != -1 &&
        (meshes[minI].mirror_ref.x > 0 || meshes[minI].mirror_ref.y > 0 || meshes[minI].mirror_ref.z > 0)) {
        mirror_existed = 1;
                        

    }
    if (mirror_existed == 1 && max_recursion > 0 && minI != -1 && minI != 90000) {
        
        Ray mirror_ray;
        Vector_3d new_camera_formirror = (this->_origin) + (P * (-1));
        new_camera_formirror = new_camera_formirror.normalize();
        float degree = N * new_camera_formirror;
        mirror_ray._direction = (new_camera_formirror * (-1)) + (N * (2 * degree));
        mirror_ray._direction=mirror_ray._direction.normalize();
        mirror_ray._origin = P + (mirror_ray._direction * (epsilon));
        for (int i = 0; i < num_spheres; i++) {
            distance_value = s[i].intersect(mirror_ray);
            if (distance_value < minT_formirror && distance_value > 0) {
                material_id_selector_formirror = 1;
                minT_formirror = distance_value;
                minI_formirror = i;
            }
        }
        for (int i = 0; i < num_triangles; i++) {
            distance_value = (triangles[i]).intersect(mirror_ray);
            if (distance_value < minT_formirror && distance_value > 0) {
                material_id_selector_formirror = 2;
                minT_formirror = distance_value;
                minI_formirror = i;
            }
        }
        for (int i = 0; i < num_meshes; i++) {
            distance_value = (meshes[i]).intersect(mirror_ray);
            if (distance_value < minT_formirror && distance_value > 0) {
                material_id_selector_formirror = 3;
                minT_formirror = distance_value;
                minI_formirror = i;
            }
        }
        if(minI_formirror != -1){
            if ( material_id_selector == 1 ){
                mirror_ray.computecolour(meshes, num_meshes, triangles, num_triangles, s, num_spheres, light, num_light,(max_recursion - 1), epsilon);
                Vector_3d new_vector_value;
                new_vector_value.x = mirror_ray.color.x * s[minI].mirror_ref.x;
                new_vector_value.y = mirror_ray.color.y * s[minI].mirror_ref.y;
                new_vector_value.z = mirror_ray.color.z * s[minI].mirror_ref.z;
                this->color = this->color + new_vector_value;
            }
            if ( material_id_selector == 2 ) {
                mirror_ray.computecolour(meshes, num_meshes, triangles, num_triangles, s, num_spheres, light, num_light, (max_recursion - 1),epsilon);
                Vector_3d new_vector_value;
                new_vector_value.x = mirror_ray.color.x * triangles[minI].mirror_ref.x;
                new_vector_value.y = mirror_ray.color.y * triangles[minI].mirror_ref.y;
                new_vector_value.z = mirror_ray.color.z * triangles[minI].mirror_ref.z;
                this->color = this->color + new_vector_value;
            }
        
            if ( material_id_selector == 3 ) 
            {   
                mirror_ray.computecolour(meshes, num_meshes, triangles, num_triangles, s, num_spheres, light, num_light, (max_recursion - 1), epsilon);
                Vector_3d new_vector_value;
                new_vector_value.x = mirror_ray.color.x * meshes[minI].mirror_ref.x;
                new_vector_value.y = mirror_ray.color.y * meshes[minI].mirror_ref.y;
                new_vector_value.z = mirror_ray.color.z * meshes[minI].mirror_ref.z;
                this->color = this->color + new_vector_value;
            }
        }

    }

}
int main(int argc, char* argv[])
{   
    //camerarr = kameralari tutan array
    //light_array = isiklari tutan array
    //image2 = olusturdugumuz gorsel
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    std::vector <parser::Camera> cameras = scene.cameras;

    Vector_3d ambient_light(scene.ambient_light.x, scene.ambient_light.y, scene.ambient_light.z);


    // for cameras
    int i = 0;
    Camera cameraarr[cameras.size()];
    float epsilon=scene.shadow_ray_epsilon;
    int max_rec=scene.max_recursion_depth;
    int num_camera = cameras.size();
    while(i < num_camera) {// HATA BURADAAAAAAAAAAAA
        //Vector_3d normalizedGaze = (cameraarr[i].gaze);
        cameraarr[i]._position = Vector_3d(cameras[i].position.x, cameras[i].position.y, cameras[i].position.z);
        cameraarr[i].gaze = Vector_3d(cameras[i].gaze.x, cameras[i].gaze.y, cameras[i].gaze.z).normalize();
        cameraarr[i]._height = cameras[i].image_height;
        cameraarr[i]._width = cameras[i].image_width;
        cameraarr[i]._distance = cameras[i].near_distance;
        cameraarr[i].nearPlane.x = cameras[i].near_plane.x;
        cameraarr[i].nearPlane.y = cameras[i].near_plane.y;
        cameraarr[i].nearPlane.w = cameras[i].near_plane.z;
        cameraarr[i].nearPlane.z = cameras[i].near_plane.w;
        cameraarr[i].image_name = cameras[i].image_name;
        cameraarr[i].up_3d = Vector_3d(cameras[i].up.x, cameras[i].up.y, cameras[i].up.z).normalize();
        cameraarr[i].up_w = ((cameraarr[i].gaze)* (-1)).normalize();
        cameraarr[i].up_u = (cameraarr[i].up_3d.cross(cameraarr[i].up_w)).normalize();
        i++;
    }

    //meshes
    std::vector <parser::Mesh> meshes = scene.meshes;
    Mesh mesh_array[meshes.size()];
        int num_mesh = meshes.size();

    for (int i = 0; i < meshes.size(); i++) {
        
        mesh_array[i].ambient_ref = Vector_3d((scene.materials[(meshes[i].material_id) - 1].ambient.x),(scene.materials[(meshes[i].material_id) - 1].ambient.y),(scene.materials[(meshes[i].material_id) - 1].ambient.z));
        mesh_array[i].spec_ref = Vector_3d((scene.materials[(meshes[i].material_id) - 1].specular.x),(scene.materials[(meshes[i].material_id) - 1].specular.y),(scene.materials[(meshes[i].material_id) - 1].specular.z));
        mesh_array[i].diffuse_ref = Vector_3d((scene.materials[(meshes[i].material_id) - 1].diffuse.x),(scene.materials[(meshes[i].material_id) - 1].diffuse.y),(scene.materials[(meshes[i].material_id) - 1].diffuse.z));
        mesh_array[i].mr_exp = scene.materials[(meshes[i].material_id) - 1].phong_exponent;
        mesh_array[i].mirror_ref = Vector_3d((scene.materials[(meshes[i].material_id) - 1].mirror.x),(scene.materials[(meshes[i].material_id) - 1].mirror.y),(scene.materials[(meshes[i].material_id) - 1].mirror.z));
        
        for (int j = 0; j < (meshes[i].faces.size()); j++) {
            Triangle tri;
            tri.first = Vector_3d((scene.vertex_data[(meshes[i].faces[j].v0_id) - 1].x),(scene.vertex_data[(meshes[i].faces[j].v0_id) - 1].y),(scene.vertex_data[(meshes[i].faces[j].v0_id) - 1].z));
            tri.second = Vector_3d((scene.vertex_data[(meshes[i].faces[j].v1_id) - 1].x),(scene.vertex_data[(meshes[i].faces[j].v1_id) - 1].y),(scene.vertex_data[(meshes[i].faces[j].v1_id) - 1].z));
            tri.third = Vector_3d((scene.vertex_data[(meshes[i].faces[j].v2_id) - 1].x),(scene.vertex_data[(meshes[i].faces[j].v2_id) - 1].y),(scene.vertex_data[(meshes[i].faces[j].v2_id) - 1].z));
            tri.ambient_ref = Vector_3d((scene.materials[(meshes[i].material_id) - 1].ambient.x),(scene.materials[(meshes[i].material_id) - 1].ambient.y),(scene.materials[(meshes[i].material_id) - 1].ambient.z));
            tri.spec_ref = Vector_3d((scene.materials[(meshes[i].material_id) - 1].specular.x),(scene.materials[(meshes[i].material_id) - 1].specular.y),(scene.materials[(meshes[i].material_id) - 1].specular.z));
            tri.diffuse_ref = Vector_3d((scene.materials[(meshes[i].material_id) - 1].diffuse.x),(scene.materials[(meshes[i].material_id) - 1].diffuse.y),(scene.materials[(meshes[i].material_id) - 1].diffuse.z));
            tri.mr_exp = scene.materials[(meshes[i].material_id) - 1].phong_exponent;
            tri.normal = ((tri.second) + (tri.first) * (-1)).cross((tri.third) + (tri.second) * (-1));
            tri.ambient_light = ambient_light;
            mesh_array[i].triangles.push_back(tri);

        }
        mesh_array[i].ambient_light = ambient_light;

    }
//spheres
    std::vector <parser::Sphere> spheres = scene.spheres;
    int num_sphere = spheres.size();

    Sphere sphere_array[spheres.size()];
    for (int i = 0; i < num_sphere; i++) {
        sphere_array[i].radius = spheres[i].radius;
        sphere_array[i].center = Vector_3d((scene.vertex_data[(spheres[i].center_vertex_id) - 1].x),(scene.vertex_data[(spheres[i].center_vertex_id) - 1].y),(scene.vertex_data[(spheres[i].center_vertex_id) - 1].z));
        sphere_array[i].ambient_ref = Vector_3d((scene.materials[(spheres[i].material_id) - 1].ambient.x),(scene.materials[(spheres[i].material_id) - 1].ambient.y),(scene.materials[(spheres[i].material_id) - 1].ambient.z));
        sphere_array[i].spec_ref = Vector_3d((scene.materials[(spheres[i].material_id) - 1].specular.x),(scene.materials[(spheres[i].material_id) - 1].specular.y),(scene.materials[(spheres[i].material_id) - 1].specular.z));
        sphere_array[i].diffuse_ref = Vector_3d((scene.materials[(spheres[i].material_id) - 1].diffuse.x),(scene.materials[(spheres[i].material_id) - 1].diffuse.y),(scene.materials[(spheres[i].material_id) - 1].diffuse.z));
        sphere_array[i].mr_exp = scene.materials[(spheres[i].material_id) - 1].phong_exponent;
        sphere_array[i].ambient_light = ambient_light;
        sphere_array[i].mirror_ref = Vector_3d((scene.materials[(spheres[i].material_id) - 1].mirror.x),(scene.materials[(spheres[i].material_id) - 1].mirror.y),(scene.materials[(spheres[i].material_id) - 1].mirror.z));
    }

    // triangle
    std::vector <parser::Triangle> triangles = scene.triangles;
    int num_triangle = triangles.size();
    Triangle triangle_array[triangles.size()];
    for (int i = 0; i < num_triangle; i++) {
        triangle_array[i].first = Vector_3d((scene.vertex_data[(scene.triangles[i].indices.v0_id) - 1].x),(scene.vertex_data[(scene.triangles[i].indices.v0_id) - 1].y),(scene.vertex_data[(scene.triangles[i].indices.v0_id) - 1].z));
        triangle_array[i].second = Vector_3d((scene.vertex_data[(scene.triangles[i].indices.v1_id) - 1].x),(scene.vertex_data[(scene.triangles[i].indices.v1_id) - 1].y),(scene.vertex_data[(scene.triangles[i].indices.v1_id) - 1].z));
        triangle_array[i].third = Vector_3d((scene.vertex_data[(scene.triangles[i].indices.v2_id) - 1].x),(scene.vertex_data[(scene.triangles[i].indices.v2_id) - 1].y),(scene.vertex_data[(scene.triangles[i].indices.v2_id) - 1].z));
        triangle_array[i].ambient_ref = Vector_3d((scene.materials[(triangles[i].material_id) - 1].ambient.x),(scene.materials[(triangles[i].material_id) - 1].ambient.y),(scene.materials[(triangles[i].material_id) - 1].ambient.z));
        triangle_array[i].spec_ref = Vector_3d((scene.materials[(triangles[i].material_id) - 1].specular.x),(scene.materials[(triangles[i].material_id) - 1].specular.y),(scene.materials[(triangles[i].material_id) - 1].specular.z));
        triangle_array[i].diffuse_ref = Vector_3d((scene.materials[(triangles[i].material_id) - 1].diffuse.x),(scene.materials[(triangles[i].material_id) - 1].diffuse.y),(scene.materials[(triangles[i].material_id) - 1].diffuse.z));
        triangle_array[i].mr_exp = scene.materials[(triangles[i].material_id) - 1].phong_exponent;
        triangle_array[i].normal = ((triangle_array[i].second) + (triangle_array[i].first) * (-1)).cross((triangle_array[i].third) + (triangle_array[i].second) * (-1));
        triangle_array[i].normal=triangle_array[i].normal.normalize();
        triangle_array[i].ambient_light = ambient_light;
        triangle_array[i].mirror_ref = Vector_3d((scene.materials[(triangles[i].material_id) - 1].mirror.x),(scene.materials[(triangles[i].material_id) - 1].mirror.y),(scene.materials[(triangles[i].material_id) - 1].mirror.z));
    }

    //for lights
    Vector_3d background_color(scene.background_color.x,scene.background_color.y,scene.background_color.z);
    std::vector <parser::PointLight> lights = scene.point_lights;

    int num_light = lights.size();
    Light light_array[num_light];
    for (int i = 0; i < num_light; i++) {
        light_array[i].position = Vector_3d(lights[i].position.x, lights[i].position.y, lights[i].position.z);
        light_array[i]._instensity = Vector_3d(lights[i].intensity.x, lights[i].intensity.y, lights[i].intensity.z);
    }
    
    for (int x = 0; x < num_camera; x++) {
        unsigned char *image2 = new unsigned char[cameraarr[x]._width * cameraarr[x]._height*3];
        int _variable = 0;
        for (int j = 0; j < cameraarr[x]._height; j++) {
            for (int i = 0; i < cameraarr[x]._width; i++) {
                Ray myray = cameraarr[x].generateRay(i, j,background_color);
                myray.computecolour(mesh_array,num_mesh,triangle_array,num_triangle,sphere_array,num_sphere,light_array,num_light,max_rec,epsilon);
                
                if (myray.color.x > 255) myray.color.x = 255;
                else if (myray.color.x < 0) myray.color.x = 0;

                image2[_variable] = (unsigned char) ((int)(myray.color.x)+0.5);

                if (myray.color.y > 255) myray.color.y = 255;
                else if (myray.color.y < 0) myray.color.y = 0;

                image2[_variable+1] = (unsigned char) ((int)(myray.color.y)+0.5);

                if (myray.color.z > 255) myray.color.z = 255;
                else if (myray.color.z < 0) myray.color.z = 0;

                image2[_variable+2] = (unsigned char) ((int)(myray.color.z)+0.5);
                _variable = _variable+3;
            }
        }
        write_ppm((cameraarr[x].image_name).c_str(), image2, cameraarr[x]._width, cameraarr[x]._height);
    }

}

