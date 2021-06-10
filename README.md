# BSPTree

This project has code from one of my last graphics projects. It divides all geometry in a scene and builds a BSP (Binary space partitioning) tree.

In BSPTree.cpp/h there is code for actually building the tree recursively using different splitting planes. 

All the other files are different types of bouding volumes to change to be more accruate or faster intersection tests.

Below is an example of what an object look like divvided into the tree.

![SubDividedObject](https://user-images.githubusercontent.com/60011821/121587540-ed9a7480-ca02-11eb-9202-5b70a1a23182.jpg)


Below are examples of what an entire scene looks like with different bounding volumes

![SphereBoundingVolume](https://user-images.githubusercontent.com/60011821/121587919-4f5ade80-ca03-11eb-8742-113e84031ca2.jpg)

![AABB](https://user-images.githubusercontent.com/60011821/121588060-74e7e800-ca03-11eb-8028-bfa78789468f.jpg)





