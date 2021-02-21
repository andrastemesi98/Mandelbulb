#include "MyApp.h"
#include "GLUtils.hpp"

#include <math.h>

CMyApp::CMyApp(void)
{
	m_vaoID = 0;
	m_vboID = 0;
	m_ibID = 0;

	m_programID = 0;

	m_camera.SetView(glm::vec3(1.5, 1.5, 0), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
}


CMyApp::~CMyApp(void)
{
}


bool CMyApp::Init()
{
	glClearColor(0.125f, 0.25f, 0.5f, 1.0f);

	glEnable(GL_CULL_FACE); 
	glEnable(GL_DEPTH_TEST);
	glCullFace(GL_BACK); 

	Vertex vert[4];
	vert[0].p = glm::vec2(-1, -1);
	vert[1].p = glm::vec2(-1, 1);
	vert[2].p = glm::vec2(1, -1);
	vert[3].p = glm::vec2(1, 1);


    GLushort indices[6];
	indices[0] = 0;
	indices[1] = 2;
	indices[2] = 3;
	indices[3] = 1;
	indices[4] = 0;
	indices[5] = 3;

	glGenVertexArrays(1, &m_vaoID);
	glBindVertexArray(m_vaoID);
	
	glGenBuffers(1, &m_vboID); 
	glBindBuffer(GL_ARRAY_BUFFER, m_vboID); 
	glBufferData( GL_ARRAY_BUFFER,	
				  sizeof(vert),		
				  vert,	
				  GL_STATIC_DRAW);	
	

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(
		0,				
		3,				
		GL_FLOAT,		
		GL_FALSE,		
		sizeof(Vertex),	
		0				
	); 

	glEnableVertexAttribArray(1); 
	glVertexAttribPointer(
		1,
		3, 
		GL_FLOAT,
		GL_FALSE,
		sizeof(Vertex),
		(void*)(sizeof(glm::vec3)) );

	glGenBuffers(1, &m_ibID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	
	GLuint vs_ID = loadShader(GL_VERTEX_SHADER,		"myVert.vert");
	GLuint fs_ID = loadShader(GL_FRAGMENT_SHADER,	"myFrag.frag");

	m_programID = glCreateProgram();

	glAttachShader(m_programID, vs_ID);
	glAttachShader(m_programID, fs_ID);

	glBindAttribLocation(	m_programID,	
							0,				
							"vs_in_pos");	
	glBindAttribLocation( m_programID, 1, "vs_in_col");

	glLinkProgram(m_programID);

	GLint infoLogLength = 0, result = 0;

	glGetProgramiv(m_programID, GL_LINK_STATUS, &result);
	glGetProgramiv(m_programID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if ( GL_FALSE == result )
	{
		std::vector<char> ProgramErrorMessage( infoLogLength );
		glGetProgramInfoLog(m_programID, infoLogLength, NULL, &ProgramErrorMessage[0]);
		fprintf(stdout, "%s\n", &ProgramErrorMessage[0]);
		
		char* aSzoveg = new char[ProgramErrorMessage.size()];
		memcpy( aSzoveg, &ProgramErrorMessage[0], ProgramErrorMessage.size());

		std::cout << "[app.Init()] Sáder Huba panasza: " << aSzoveg << std::endl;

		delete aSzoveg;
	}

	glDeleteShader( vs_ID );
	glDeleteShader( fs_ID );

	//
	// egyéb inicializálás
	//

	// vetítési mátrix létrehozása
	m_matProj = glm::perspective( 45.0f, 640/480.0f, 1.0f, 1000.0f );

	// shader-beli transzformációs mátrixok címének lekérdezése
	m_loc_mvp = glGetUniformLocation(m_programID, "MVP");
	m_loc_win = glGetUniformLocation(m_programID, "mywin");
	m_loc_eye = glGetUniformLocation(m_programID, "cameye");
	m_loc_at = glGetUniformLocation(m_programID, "camat");
	m_loc_up = glGetUniformLocation(m_programID, "camup");
	m_loc_time = glGetUniformLocation(m_programID, "mytime");

	m_camera.SetProj(45.0f, 640.0f / 480.0f, 0.01f, 1000.0f);

	return true;
}

void CMyApp::Clean()
{
	glDeleteBuffers(1, &m_vboID);
	glDeleteBuffers(1, &m_ibID);
	glDeleteVertexArrays(1, &m_vaoID);

	glDeleteProgram( m_programID );
}

void CMyApp::Update()
{
	static Uint32 last_time = SDL_GetTicks();
	float delta_time = (SDL_GetTicks() - last_time) / 1000.0f;

	m_camera.Update(delta_time);

	last_time = SDL_GetTicks();
}


void CMyApp::Render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram( m_programID );


	m_matWorld = glm::mat4(1.0f);

	glm::mat4 mvp = m_matProj * m_matView * m_matWorld;

	glUniformMatrix4fv( m_loc_mvp,
						1,			
						GL_FALSE,	
						&(mvp[0][0]) );

	glm::vec2 seged = glm::vec2(width, height);
	seged = glm::normalize(seged);

	glUniform2f(m_loc_win, seged.x, seged.y);

	glUniform3f(m_loc_eye, m_camera.GetEye().x, m_camera.GetEye().y, m_camera.GetEye().z);
	glUniform3f(m_loc_at, m_camera.GetAt().x, m_camera.GetAt().y, m_camera.GetAt().z);
	glUniform3f(m_loc_up, m_camera.GetUp().x, m_camera.GetUp().y, m_camera.GetUp().z);

	glUniform1f(m_loc_time, SDL_GetTicks());

	
	glBindVertexArray(m_vaoID);

	glDrawElements(	GL_TRIANGLES,		
					6,		
					GL_UNSIGNED_SHORT,	
					0);					

	glBindVertexArray(0);

	glUseProgram( 0 );
}

void CMyApp::KeyboardDown(SDL_KeyboardEvent& key)
{
	m_camera.KeyboardDown(key);
}

void CMyApp::KeyboardUp(SDL_KeyboardEvent& key)
{
	m_camera.KeyboardUp(key);
}

void CMyApp::MouseMove(SDL_MouseMotionEvent& mouse)
{
	m_camera.MouseMove(mouse);
}

void CMyApp::MouseDown(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseUp(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseWheel(SDL_MouseWheelEvent& wheel)
{
}

void CMyApp::Resize(int _w, int _h)
{
	glViewport(0, 0, _w, _h);

	m_matProj = glm::perspective(  45.0f,		
									_w/(float)_h,	
									0.01f,			
									100.0f);

	width = _w;
	height = _h;

	m_camera.Resize(_w, _h);
}