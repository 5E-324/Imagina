#ifndef _GLUTIL_H_
#define _GLUTIL_H_

// PLATFORM DEPENDENT
extern PFNGLCREATEPROGRAMPROC glCreateProgram;
extern PFNGLCREATESHADERPROC glCreateShader;
extern PFNGLSHADERSOURCEPROC glShaderSource;
extern PFNGLCOMPILESHADERPROC glCompileShader;
extern PFNGLGETSHADERIVPROC glGetShaderiv;
extern PFNGLGETPROGRAMIVPROC glGetProgramiv;
extern PFNGLDELETESHADERPROC glDeleteShader;
extern PFNGLATTACHSHADERPROC glAttachShader;
extern PFNGLLINKPROGRAMPROC glLinkProgram;
extern PFNGLUSEPROGRAMPROC	glUseProgram;

extern PFNGLGENVERTEXARRAYSPROC			glGenVertexArrays;
extern PFNGLBINDVERTEXARRAYPROC			glBindVertexArray;
extern PFNGLGENBUFFERSPROC				glGenBuffers;
extern PFNGLBINDBUFFERPROC				glBindBuffer;
extern PFNGLBUFFERDATAPROC				glBufferData;
extern PFNGLBUFFERSUBDATAPROC			glBufferSubData;
extern PFNGLGETATTRIBLOCATIONPROC		glGetAttribLocation;
extern PFNGLENABLEVERTEXATTRIBARRAYPROC	glEnableVertexAttribArray;
extern PFNGLVERTEXATTRIBPOINTERPROC		glVertexAttribPointer;
extern PFNGLACTIVETEXTUREPROC			glActiveTexture;
extern PFNGLGETUNIFORMLOCATIONPROC		glGetUniformLocation;
extern PFNGLUNIFORM1IPROC				glUniform1i;
extern PFNGLUNIFORM1FPROC				glUniform1f;

extern PFNWGLSWAPINTERVALEXTPROC		wglSwapIntervalEXT;

extern PFNGLGETSHADERINFOLOGPROC		glGetShaderInfoLog;
extern PFNGLGETPROGRAMINFOLOGPROC		glGetProgramInfoLog;

void InitGLExtensions();

GLuint LoadShaders(const char *VSSource, const char *FSSource);

#endif