#pragma once

#include <QVector3D>
#include <QVector4D>

namespace Waves {
	class Vertex
	{
	public:
		// These members are public since we don't need to do do anything special with regards to accessors/mutators
		QVector3D position{ 0.0,0.0,0.0 };
		QVector3D normal{ 0.0,0.0,1.0 };
		float shininess{ 0.8f };
		QVector4D specular{ 1.0,1.0,1.0,1.0 };

		// Constructors
		Vertex() = default;
		Vertex(float x, float y) : position({ x, y, 0.0f }) {}
		Vertex(const QVector3D& position) : position(position) {}
		Vertex(const QVector3D& position, const QVector3D& normal) : position(position), normal(normal) {}
		Vertex(const QVector3D& position, const QVector3D& normal, const float& shininess, const QVector4D& specular) : position(position), normal(normal), shininess(shininess), specular(specular) {}

		// OpenGL Helpers
		static const int PositionTupleSize = 3;
		static const int NormalTupleSize = 3;
		static const int ShininessTupleSize = 1;
		static const int SpecularTupleSize = 4;
		static constexpr int positionOffset() { return offsetof(Vertex, position); }
		static constexpr int normalOffset() { return offsetof(Vertex, normal); }
		static constexpr int shininessOffset() { return offsetof(Vertex, shininess); }
		static constexpr int specularOffset() { return offsetof(Vertex, specular); }
		static constexpr int stride() { return sizeof(Vertex); }

	};
}