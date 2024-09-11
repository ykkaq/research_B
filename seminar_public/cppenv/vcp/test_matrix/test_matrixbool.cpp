// VCP Library
// http ://verified.computation.jp
//   
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
// Copyright(c) 2017, Kouta Sekine <k.sekine@computation.jp>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met :
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and / or other materials provided with the distribution.
// * Neither the name of the Kouta Sekine nor the names of its contributors
//   may be used to endorse or promote products derived from this software
//   without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED.IN NO EVENT SHALL KOUTA SEKINE BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <iostream>

#include <omp.h>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/imats.hpp>
#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>


int main(void) {
	 vcp::matrix< bool > A, B, C;
	 //vcp::mbool A, B, C;

	A.allfalse(1, 10);
	A(0, 1) = true;
	B.allfalse(1, 10);
	B(0, 1) = true;
	B(0, 0) = true;

	std::cout << A << std::endl;
	std::cout << B << std::endl;

	C = A && B;
	std::cout << C << std::endl;

	C = !((A || B) || B);
	std::cout << C << std::endl;
	std::cout << !C(0, 1) << std::endl;

	C.resize(5, 11);
	std::cout << C << std::endl;

	C.clear();

	vcp::matrix< kv::interval< double >, vcp::pidblas > AA, BB, CC;
	AA.rand(10);
	BB.rand(10);
	CC.rand(10);
	std::cout << ((AA != BB) || (AA != CC)) << std::endl;;

	C(0, 0).flip();
	C = AA != BB;

	std::cout << all((C || (AA != BB)) && (AA == CC)) << std::endl;
	std::cout << any((C || (AA != BB)) && (AA == CC)) << std::endl;
	std::cout << none((C || (AA != BB)) && (AA == CC)) << std::endl;

	if (none((C || (AA != BB)) && (AA == CC))) {
		std::cout << "nya---" << std::endl;
	}
}