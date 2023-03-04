#include "Includes.h"
#include <iostream>

PFNGLCREATEPROGRAMPROC	glCreateProgram;
PFNGLCREATESHADERPROC	glCreateShader;
PFNGLSHADERSOURCEPROC	glShaderSource;
PFNGLCOMPILESHADERPROC	glCompileShader;
PFNGLGETSHADERIVPROC	glGetShaderiv;
PFNGLGETPROGRAMIVPROC	glGetProgramiv;
PFNGLDELETESHADERPROC	glDeleteShader;
PFNGLATTACHSHADERPROC	glAttachShader;
PFNGLLINKPROGRAMPROC	glLinkProgram;
PFNGLUSEPROGRAMPROC		glUseProgram;

PFNGLGENVERTEXARRAYSPROC			glGenVertexArrays;
PFNGLBINDVERTEXARRAYPROC			glBindVertexArray;
PFNGLGENBUFFERSPROC					glGenBuffers;
PFNGLBINDBUFFERPROC					glBindBuffer;
PFNGLBUFFERDATAPROC					glBufferData;
PFNGLBUFFERSUBDATAPROC				glBufferSubData;
PFNGLGETATTRIBLOCATIONPROC			glGetAttribLocation;
PFNGLENABLEVERTEXATTRIBARRAYPROC	glEnableVertexAttribArray;
PFNGLVERTEXATTRIBPOINTERPROC		glVertexAttribPointer;
PFNGLACTIVETEXTUREPROC				glActiveTexture;
PFNGLGETUNIFORMLOCATIONPROC			glGetUniformLocation;
PFNGLUNIFORM1IPROC					glUniform1i;
PFNGLUNIFORM1FPROC					glUniform1f;

PFNWGLSWAPINTERVALEXTPROC			wglSwapIntervalEXT;

PFNGLGETSHADERINFOLOGPROC			glGetShaderInfoLog;
PFNGLGETPROGRAMINFOLOGPROC			glGetProgramInfoLog;

void InitGLExtensions() {
	glCreateProgram	= (PFNGLCREATEPROGRAMPROC)	wglGetProcAddress("glCreateProgram");
	glCreateShader	= (PFNGLCREATESHADERPROC)	wglGetProcAddress("glCreateShader");
	glShaderSource	= (PFNGLSHADERSOURCEPROC)	wglGetProcAddress("glShaderSource");
	glCompileShader	= (PFNGLCOMPILESHADERPROC)	wglGetProcAddress("glCompileShader");
	glGetShaderiv	= (PFNGLGETSHADERIVPROC)	wglGetProcAddress("glGetShaderiv");
	glGetProgramiv	= (PFNGLGETPROGRAMIVPROC)	wglGetProcAddress("glGetProgramiv");
	glDeleteShader	= (PFNGLDELETESHADERPROC)	wglGetProcAddress("glDeleteShader");
	glAttachShader	= (PFNGLATTACHSHADERPROC)	wglGetProcAddress("glAttachShader");
	glLinkProgram	= (PFNGLLINKPROGRAMPROC)	wglGetProcAddress("glLinkProgram");
	glUseProgram	= (PFNGLUSEPROGRAMPROC)		wglGetProcAddress("glUseProgram");
	
	glGenVertexArrays			= (PFNGLGENVERTEXARRAYSPROC)			wglGetProcAddress("glGenVertexArrays");
	glBindVertexArray			= (PFNGLBINDVERTEXARRAYPROC)			wglGetProcAddress("glBindVertexArray");
	glGenBuffers				= (PFNGLGENBUFFERSPROC)					wglGetProcAddress("glGenBuffers");
	glBindBuffer				= (PFNGLBINDBUFFERPROC)					wglGetProcAddress("glBindBuffer");
	glBufferData				= (PFNGLBUFFERDATAPROC)					wglGetProcAddress("glBufferData");
	glBufferSubData				= (PFNGLBUFFERSUBDATAPROC)				wglGetProcAddress("glBufferSubData");
	glGetAttribLocation			= (PFNGLGETATTRIBLOCATIONPROC)			wglGetProcAddress("glGetAttribLocation");
	glEnableVertexAttribArray	= (PFNGLENABLEVERTEXATTRIBARRAYPROC)	wglGetProcAddress("glEnableVertexAttribArray");
	glVertexAttribPointer		= (PFNGLVERTEXATTRIBPOINTERPROC)		wglGetProcAddress("glVertexAttribPointer");
	glActiveTexture				= (PFNGLACTIVETEXTUREPROC)				wglGetProcAddress("glActiveTexture");
	glGetUniformLocation		= (PFNGLGETUNIFORMLOCATIONPROC)			wglGetProcAddress("glGetUniformLocation");
	glUniform1i					= (PFNGLUNIFORM1IPROC)					wglGetProcAddress("glUniform1i");
	glUniform1f					= (PFNGLUNIFORM1FPROC)					wglGetProcAddress("glUniform1f");

	wglSwapIntervalEXT			= (PFNWGLSWAPINTERVALEXTPROC)			wglGetProcAddress("wglSwapIntervalEXT");

	glGetShaderInfoLog			= (PFNGLGETSHADERINFOLOGPROC)			wglGetProcAddress("glGetShaderInfoLog");
	glGetProgramInfoLog			= (PFNGLGETPROGRAMINFOLOGPROC)			wglGetProcAddress("glGetProgramInfoLog");
}

GLuint LoadShaders(const char *VSSource, const char *FSSource) {
	GLuint Program = glCreateProgram();

	GLuint VertexShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(VertexShader, 1, &VSSource, nullptr);
	glShaderSource(FragmentShader, 1, &FSSource, nullptr);
	glCompileShader(VertexShader);

	GLint Compiled;
	glGetShaderiv(VertexShader, GL_COMPILE_STATUS, &Compiled);
	if (!Compiled) {
#ifdef _DEBUG
		GLsizei Length;
		glGetShaderiv(VertexShader, GL_INFO_LOG_LENGTH, &Length);

		GLchar *Log = new GLchar[Length + 1];
		glGetShaderInfoLog(VertexShader, Length, &Length, Log);
		std::cerr << "Shader compilation failed: " << Log << std::endl;
		delete[] Log;
#endif /* DEBUG */
		glDeleteShader(VertexShader);
		return 0;
	}

	glCompileShader(FragmentShader);

	glGetShaderiv(FragmentShader, GL_COMPILE_STATUS, &Compiled);
	if (!Compiled) {
#ifdef _DEBUG
		GLsizei Length;
		glGetShaderiv(FragmentShader, GL_INFO_LOG_LENGTH, &Length);

		GLchar *Log = new GLchar[Length + 1];
		glGetShaderInfoLog(FragmentShader, Length, &Length, Log);
		std::cerr << "Shader compilation failed: " << Log << std::endl;
		delete[] Log;
#endif /* DEBUG */
		glDeleteShader(VertexShader);
		glDeleteShader(FragmentShader);
		return 0;
	}

	glAttachShader(Program, VertexShader);
	glAttachShader(Program, FragmentShader);

	glLinkProgram(Program);

	GLint Linked;
	glGetProgramiv(Program, GL_LINK_STATUS, &Linked);
	if (!Linked) {
#ifdef _DEBUG
		GLsizei Length;
		glGetProgramiv(Program, GL_INFO_LOG_LENGTH, &Length);

		GLchar *Log = new GLchar[Length + 1];
		glGetProgramInfoLog(Program, Length, &Length, Log);
		std::cerr << "Shader linking failed: " << Log << std::endl;
		delete[] Log;
#endif /* DEBUG */

		glDeleteShader(VertexShader);
		glDeleteShader(FragmentShader);

		return 0;
	}

	return Program;
}