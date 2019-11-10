"use strict";

/**
 * This class implements frequently used numerical algorithms such as the inverse power method.
 * @memberof module:Utils
 */
class Solvers {
	/**
	 * Computes the residual of Ax - λx, where x has unit norm and λ = x.Ax.
	 * @param {module:LinearAlgebra.ComplexSparseMatrix} A The complex sparse matrix whose eigen decomposition
	 * is being computed.
	 * @param {module:LinearAlgebra.ComplexDenseMatrix} x The current guess for the smallest eigenvector
	 * (corresponding to the smallest eigenvalue λ) of A.
	 * @returns {number}
	 */
	static residual(A, x) {
		x.scaleBy(new Complex(1/x.norm(),0))
		let Ax = A.timesDense(x)
		let eig = x.transpose().conjugate().timesDense(Ax)
		let eigenVal = eig.get(0,0)
		let lambdaX = x.timesComplex(eigenVal)
		let residVec = Ax.minus(lambdaX)
		let resid = residVec.norm()
		return resid;
	}

	/**
	 * Solves Ax = λx, where λ is the smallest nonzero eigenvalue of A and x is the
	 * corresponding eigenvector. x should be initialized to a random complex dense
	 * vector (i.e., x.nCols() == 1).
	 * @param {module:LinearAlgebra.ComplexSparseMatrix} A The complex positive definite sparse matrix
	 * whose eigen decomposition needs to be computed.
	 * @returns {module:LinearAlgebra.ComplexDenseMatrix} The smallest eigenvector (corresponding to the
	 * smallest eigenvalue λ) of A.
	 */
	static solveInversePowerMethod(A) {
		let x = ComplexDenseMatrix.random(A.nRows(), 1)
		x.scaleBy(new Complex(1/x.norm()))
		const maxIters = 1000
		let iters = 0
		while(this.residual(A, x) > 1e-10 && iters < maxIters) {
			iters++
			x = A.lu().solveSquare(x)
			let xBar = x.sum().overReal(x.nRows())
			let xBarVec = ComplexDenseMatrix.constant(xBar, x.nRows(), 1)
			x.decrementBy(xBarVec)
			x.scaleBy(new Complex(1/x.norm(2)))
		}
		this.residual(A, x)
		console.log(x.get(0,0), x.get(1,0), x.get(2,0))
		let test = A.timesDense(x)
		test.scaleBy(new Complex(1/test.norm()))
		console.log(test.get(0,0), test.get(1,0), test.get(2,0))
		return x
	}

	/**
	 * Inverts a 2x2 matrix.
	 * @param {module:LinearAlgebra.DenseMatrix} m The matrix to be inverted.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	static invert2x2(m) {
		let m00 = m.get(0, 0);
		let m01 = m.get(0, 1);
		let m10 = m.get(1, 0);
		let m11 = m.get(1, 1);

		let det = m00 * m11 - m01 * m10;
		m.set(m11, 0, 0);
		m.set(m00, 1, 1);
		m.set(-m01, 0, 1);
		m.set(-m10, 1, 0);
		m.scaleBy(1.0 / det);

		return m;
	}
}