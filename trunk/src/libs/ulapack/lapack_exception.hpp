/* uLapack exception classes.
 *
 * uLapack functions can throw two types of exception: NumericalError and
 * LogicalError. The first indicates a computational failure and the latter
 * indicates a logical flaw in the program (such as incorrect function arguments).
 * Both exceptions are derived from LapackError, which in turn is derived from
 * std::exception, so all uLapack exceptions can be caught using these common
 * base types.
 *
 * Note:
 * (1) All LAPACK functions update an integer "info" value:
 *	info == 0 : successful termination 
 *	info < 0  : illegal value of one or more arguments -- no computation performed
 *	info > 0  : failure in the course of computation 
 * If info < 0, the value of -info indicates the argument in the function argument
 * list that is invalid. If info > 0, consult the lapack documentation for the error
 * code.
 *
 * (2) All exceptions derived from LapackError have an "info" value that can be
 * obtained by calling the get_info() member. This value is 0 for non-clapack errors.
 *
 * Tim Bailey 2005.
 */

#ifndef ULAPACK_EXCEPTIONS_HPP_
#define ULAPACK_EXCEPTIONS_HPP_

#include <exception>

namespace ulapack {

	class LapackError : public std::exception
	{
		int info_;
		const char *message_;
	protected:
		LapackError(const char *message, int info) : 
			 info_(info), message_(message) {}
	public:
		int get_info() const { return info_; }
		const char* what() const throw() { return message_; }
	};

	class NumericalError : public LapackError
	{
	public:
		NumericalError(const char *message, int info=0) :
		  LapackError(message, info) {}
	};

	class LogicalError : public LapackError
	{
	public:
		LogicalError(const char *message, int info=0) :
		  LapackError(message, info) {}
	};

} // namespace ulapack

#endif
