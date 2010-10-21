		#-Dcimg_use_zlib -Dcimg_use_magick -Dcimg_use_fftw3 -Dcimg_use_lapack \
		#-Dcimg_use_openmp -Dcimg_use_png -Dcimg_use_jpeg -Dcimg_use_tiff -Dcimg_use_ffmpeg \
		#-fopenmp -Dcimg_use_openmp \

CXXFLAGS =	-O3 -Wall -Wno-deprecated \
		-Dcimg_use_vt100 -Dcimg_use_xshm -Dcimg_use_xrandr \
		-Dcimg_use_png \
		-I/usr/X11R6/include 

OBJS =		voronoi.o PolyLine_Simplification.o

LIBS = -L/usr/X11R6/lib -lqhull -lpthread -lX11 -lXext -lXrandr -lpng -g

TARGET =	voronoi

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
