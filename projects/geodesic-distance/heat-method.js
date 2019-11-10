"use strict";

class HeatMethod {
	/**
	 * This class implements the {@link http://cs.cmu.edu/~kmcrane/Projects/HeatMethod/ heat method} to compute geodesic distance
	 * on a surface mesh.
	 * @constructor module:Projects.HeatMethod
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} F The mean curvature flow operator built on the input mesh.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		const vertices = geometry.mesh.vertices
		this.vertexIndex = indexElements(vertices);

		this.A = this.geometry.laplaceMatrix(this.vertexIndex ) 

		let h = 0
		// calculate mean spacing
		for(let e of this.geometry.mesh.edges) {
			h += this.geometry.length(e)
		}
		h = h / this.geometry.mesh.edges.length
		h = h * h
		// flow op
		let L = this.A.timesReal(1.0)
		L.scaleBy(h)
		const M = this.geometry.massMatrix(this.vertexIndex)
		let op = M.plus(L)
		this.F = op
	}

	/**
	 * Computes the vector field X = -∇u / |∇u|.
	 * @private
	 * @method module:Projects.HeatMethod#computeVectorField
	 * @param {module:LinearAlgebra.DenseMatrix} u A dense vector (i.e., u.nCols() == 1) representing the
	 * heat that is allowed to diffuse on the input mesh for a brief period of time.
	 * @returns {Object} A dictionary mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 */
	computeVectorField(u) {
		let X = {}
		for(let f of this.geometry.mesh.faces){
			const N = this.geometry.faceNormal(f)
			let total = new Vector()
			for(let c of f.adjacentCorners()){
				const v = c.vertex
				const i = this.vertexIndex[v]
				const ei = this.geometry.vector(c.halfedge)
				let res = ei.cross(N)
				const ui = u.get(i, 0)
				res.scaleBy(ui)
				total.incrementBy(res)
			}
			total.normalize()
			X[f] = total
		}
		return X
	}

	/**
	 * Computes the integrated divergence ∇.X.
	 * @private
	 * @method module:Projects.HeatMethod#computeDivergence
	 * @param {Object} X The vector field -∇u / |∇u| represented by a dictionary
	 * mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	computeDivergence(X) {
		let gradX = DenseMatrix.zeros(this.geometry.mesh.vertices.length, 1)
		for(let v of this.geometry.mesh.vertices) {
			let div = 0
			for(let h0 of v.adjacentHalfedges()){
				const f = h0.face
				const Xj = X[f]
				let hs = [h0, h0.next.next]
				let e = hs.map(h => this.geometry.vector(h))
				e[1].scaleBy(-1)
				let cot = hs.map(h => this.geometry.cotan(h))
				let eXj = e.map(ei => Xj.dot(ei))
				let eXjCot = eXj.map((e,i) => e * cot[i])
				div += eXjCot.reduce((a,b)=> a+b)
			}
			div = div/2
			gradX.set(div, this.vertexIndex[v], 0)
		}
		return gradX
	}

	/**
	 * Shifts φ such that its minimum value is zero.
	 * @private
	 * @method module:Projects.HeatMethod#subtractMinimumDistance
	 * @param {module:LinearAlgebra.DenseMatrix} phi The (minimum 0) solution to the poisson equation Δφ = ∇.X.
	 */
	subtractMinimumDistance(phi) {
		let min = Infinity;
		for (let i = 0; i < phi.nRows(); i++) {
			min = Math.min(phi.get(i, 0), min);
		}

		for (let i = 0; i < phi.nRows(); i++) {
			phi.set(phi.get(i, 0) - min, i, 0);
		}
	}

	/**
	 * Computes the geodesic distances φ using the heat method.
	 * @method module:Projects.HeatMethod#compute
	 * @param {module:LinearAlgebra.DenseMatrix} delta A dense vector (i.e., delta.nCols() == 1) containing
	 * heat sources, i.e., u0 = δ(x).
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	compute(delta) {
		let u = this.F.chol().solvePositiveDefinite(delta)
		let X = this.computeVectorField(u)
		let divX = this.computeDivergence(X)
		let phi = this.A.chol().solvePositiveDefinite(divX)
		phi.scaleBy(-1)
		// since φ is unique up to an additive constant, it should
		// be shifted such that the smallest distance is zero
		this.subtractMinimumDistance(phi);

		return phi;
	}
}
