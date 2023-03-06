prefix = $(HOME)/win/posix/x86_64
MPIR = $(HOME)/win/src/mpir
OPENGL= $(HOME)/win/src/OpenGL-Registry

CPPFLAGS = -I$(prefix)/include -I$(MPIR) -I$(OPENGL)/api/GL
LDFLAGS = -L$(prefix)/lib -static

CXXFLAGS = -MMD -std=c++20 -mavx2 -mfma -fno-math-errno -flto -O3

LIBS = -lmpir -lmpirxx -lpng -lz -lcomctl32 -lcomdlg32 -lgdi32 -lopengl32

CXX = x86_64-w64-mingw32-g++
WINDRES = x86_64-w64-mingw32-windres

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES)) Resource.o
DEPENDS = $(patsubst %.o,%.d,$(OBJECTS))

all: Imagina.exe

clean:
	-rm -f $(OBJECTS) $(DEPENDS) Resource.res Imagina.exe

Imagina.exe: $(OBJECTS)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o Imagina.exe $(OBJECTS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

Resource.o: Resource.res
	$(WINDRES) -J res -i Resource.res -o Resource.o

Resource.res: Resource.rc
	$(WINDRES) -i Resource.rc -o Resource.res

-include $(DEPENDS)
