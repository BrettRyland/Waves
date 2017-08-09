#pragma once
///@file

#include <QVector3D>
#include <QVector4D>

namespace Waves {
	/** A class containing the vertex information that OpenGL requires.
	@note If this class is changed, then the relevant changes must be propagated to the OpenCL kernel in mesh.cl.
	*/
	class Vertex
	{
	public:
		// These members are public since we don't need to do do anything special with regards to accessors/mutators
		QVector3D position{ 0.0f, 0.0f, 0.0f }; ///< Position in R^3
		QVector3D normal{ 0.0f, 0.0f, 1.0f }; ///< Vertex normal in R^3
		float shininess{ 20.0f }; ///< Material shininess factor
		QVector4D specular{ 0.1f, 0.3f, 0.7f, 0.1f }; ///< Material specular factor

		/** @name Constructors and Destructors
		@{ */
		/// Constructor with position, normal, shininess and specular components
		Vertex(const QVector3D& position, const QVector3D& normal, const float& shininess, const QVector4D& specular) : position(position), normal(normal), shininess(shininess), specular(specular) {}
		Vertex(const QVector3D& position, const QVector3D& normal) : position(position), normal(normal) {} ///< \overload
		Vertex(const QVector3D& position) : position(position) {} ///< \overload
		Vertex(float x, float y) : position({ x, y, 0.0f }) {} ///< \overload
		Vertex() = default; ///< \overload
		///@}

		/** @name OpenGL helper values
		@{ */
		static const int PositionTupleSize = 3; ///< Number of elements in position
		static const int NormalTupleSize = 3; ///< Number of elements in normal
		static const int ShininessTupleSize = 1; ///< Number of elements in shininess
		static const int SpecularTupleSize = 4; ///< Number of elements in specular
		///@}
		/** @name OpenGL helper functions
		@{ */
		static constexpr int positionOffset() { return offsetof(Vertex, position); } ///< Offset in memory of position
		static constexpr int normalOffset() { return offsetof(Vertex, normal); } ///< Offset in memory of normal
		static constexpr int shininessOffset() { return offsetof(Vertex, shininess); } ///< Offset in memory of shininess
		static constexpr int specularOffset() { return offsetof(Vertex, specular); } ///< Offset in memory of specular
		static constexpr int stride() { return sizeof(Vertex); } ///< Size in memory of an instance of the class
		///@}
	};
}