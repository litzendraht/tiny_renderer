# Description

Rust practice project, implementing a renderer, based on the Dmitry V. Sokolov's [course](https://github.com/ssloy/tinyrenderer).

Supports several shading pipelines at the same time at a cost of slight bloat, examples below.

# Usage

Pressing `q`, `e` rotates the light, pressing `a`, `d` rotates the camera.

Binary can be launched as is to do a render of diablo with default pipeline but also very crudely supports passing of 2 arguments:

`-p`   Path to the asset folder, e.g `-p assets/african_head`.

`-s`   Choice of the shader pipeline, e.g `-s default`. All possible options:
- default
- phong
- true_normal
- darboux
- specular
- shadow
- occlusion

For now asset folder is required to contain 5 files:
- model.obj 
- normal_map_tangent.tga
- normal_map.tga
- specular_map.tga
- texture.tga (diffuse texture)

## Shading with face normals
![image](https://user-images.githubusercontent.com/17012740/211145601-f1adc7a7-fcb6-49a4-a1c6-f5b8def02b73.png)

## Phong shading
![image](https://user-images.githubusercontent.com/17012740/211145627-be524150-d663-456a-8b7f-3de8bf3d7cf6.png)

## Shading using normal map
![image](https://user-images.githubusercontent.com/17012740/211145663-1194bc4f-43a9-40df-8817-c7191922d765.png)

## Shading using specular map
Doesn't work too well for the diablo model for some reason, so here is the african head model.
![image](https://user-images.githubusercontent.com/17012740/211145711-75b00d53-98d5-4d8c-81f9-ad91de8c629d.png)

## Shading by using normals in tangent basis
Notice how tail and right hand are correctly shaded in comparison to normal map based shading. 
![image](https://user-images.githubusercontent.com/17012740/211145781-24c2a16d-9155-4c2a-b72d-98e0e1d4e504.png)

## Hard shadows + Phong shading
![image](https://user-images.githubusercontent.com/17012740/211145845-af1ac974-7168-45e1-bc2d-7445560b6dc5.png)

## Ambient occlusion
![image](https://user-images.githubusercontent.com/17012740/211145561-d249adcf-4dc5-4b5c-867a-407039a65691.png)
