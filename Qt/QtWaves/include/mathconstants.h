namespace br {
	/// Templated expressions for pi/tau, use as pi<T>/tau<T>
	template<class T> constexpr T pi = T(3.14159265358979323846264338328L);
	template<class T> constexpr T tau = pi<T>*static_cast<T>(2.0L);
}