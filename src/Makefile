EIGEN_INC_DIR=/usr/local/include/eigen3
GL_LIBS=-framework Accelerate -framework QTKit -framework CoreVideo -framework Carbon -framework GLUT -framework OpenGL -framework Cocoa -framework QuickTime -L/usr/local/lib -ljpeg
ARTK_CFLAGS= -stdlib=libc++ -O -DHAVE_NFT=1 -Wno-write-strings -Wno-deprecated-declarations -I $(ARTK_INC_DIR) -I $(EIGEN_INC_DIR)
IK_CFLAGS= -O -Wno-deprecated-declarations -I$(EIGEN_INC_DIR)

all: SuperHelix

SuperHelix: SuperHelix.cpp SuperHelix.hpp
	c++ $(IK_CFLAGS) $(GL_LIBS)  SuperHelix.cpp -o SuperHelix

