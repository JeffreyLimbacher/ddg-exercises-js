"use strict";

class HodgeDecomposition {
	/**
	 * This class computes the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf hodge decomposition} of a vector field on a surface mesh.
	 * @constructor module:Projects.HodgeDecomposition
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} edgeIndex A dictionary mapping each edge of the input mesh to a unique index.
	 * @property {module:LinearAlgebra.SparseMatrix} hodge1 The hodge star 1-form matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} hodge2 The hodge star 2-form matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} d0 The exterior derivaitve 0-form matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} d1 The exterior derivaitve 1-form matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} hodge1Inv The inverse hodge star 1-form matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} hodge2Inv The inverse hodge star 2-form matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} d0T Transpose of the exterior derivaitve 0-form matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} d1T Transpose of the exterior derivaitve 1-form matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} A The 0-form laplace matrix d0^T star1 d0 of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} B The 2-form matrix d1 star1^-1 d1^T of the input mesh.
	 */
	constructor(geometry) {
		// index vertices, edges and faces
		let vertexIndex = indexElements(geometry.mesh.vertices);
		this.edgeIndex = indexElements(geometry.mesh.edges);
		let faceIndex = indexElements(geometry.mesh.faces);

		this.hodge1 = DEC.buildHodgeStar1Form(geometry, this.edgeIndex) 
		this.hodge2 = DEC.buildHodgeStar2Form(geometry, faceIndex)
		this.d0 = DEC.buildExteriorDerivative0Form(geometry, this.edgeIndex, vertexIndex)
		this.d1 = DEC.buildExteriorDerivative1Form(geometry, faceIndex, this.edgeIndex)

		this.hodge1Inv = this.hodge1.invertDiagonal()
		this.hodge2Inv = this.hodge2.invertDiagonal()
		this.d0T = this.d0.transpose() 
		this.d1T = this.d1.transpose() 
		this.A = geometry.laplaceMatrix(vertexIndex)

		this.B = this.d1.timesSparse(this.hodge1Inv.timesSparse(this.d1T))
	}

	/**
	 * Computes the 0-form potential α by solving the system 𝛿dα = 𝛿ω.
	 * @method module:Projects.HodgeDecomposition#computeExactComponent
	 * @param {module:LinearAlgebra.DenseMatrix} omega A 1-form on the edges of the input mesh.
	 * @returns {module:LinearAlgebra.DenseMatrix} The exact component dα of ω.
	 */
	computeExactComponent(omega) {
		const omegaH = this.hodge1.timesDense(omega)
		const rhs = this.d0T.timesDense(omegaH)
		const alpha = this.A.chol().solvePositiveDefinite(rhs)
		const dAlpha = this.d0.timesDense(alpha)
		return dAlpha
	}

	/**
	 * Computes the 2-form potential β by solving the system d𝛿β = dω.
	 * @method module:Projects.HodgeDecomposition#computeCoExactComponent
	 * @param {module:LinearAlgebra.DenseMatrix} omega A 1-form on the edges of the input mesh.
	 * @returns {module:LinearAlgebra.DenseMatrix} The coexact component 𝛿β of ω.
	 */
	computeCoExactComponent(omega) {
		const rhs = this.d1.timesDense(omega)
		const betaTilde = this.B.lu().solveSquare(rhs)
		const deltaBeta = this.hodge1Inv.timesDense(this.d1T.timesDense(betaTilde))
		return deltaBeta
	}

	/**
	 * Computes the harmonic component γ = ω - dα - 𝛿β of ω.
	 * @method module:Projects.HodgeDecomposition#computeHarmonicComponent
	 * @param {module:LinearAlgebra.DenseMatrix} omega A 1-form on the edges of the input mesh.
	 * @param {module:LinearAlgebra.DenseMatrix} dAlpha The exact component of ω.
	 * @param {module:LinearAlgebra.DenseMatrix} deltaBeta The coexact component of ω.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	computeHarmonicComponent(omega, dAlpha, deltaBeta) {
		const gamma = omega.minus(dAlpha).minus(deltaBeta)
		return gamma
	}
}
