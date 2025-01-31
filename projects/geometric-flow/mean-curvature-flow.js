"use strict";

class MeanCurvatureFlow {
	/**
	 * This class performs {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf mean curvature flow} on a surface mesh.
	 * @constructor module:Projects.MeanCurvatureFlow
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the mean curvature flow operator.
	 * @private
	 * @method module:Projects.MeanCurvatureFlow#buildFlowOperator
	 * @param {module:LinearAlgebra.SparseMatrix} M The mass matrix of the input mesh.
	 * @param {number} h The timestep.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	buildFlowOperator(M, h) {
		let L = this.geometry.laplaceMatrix(this.vertexIndex)
		L.scaleBy(h)
		let op = M.plus(L)
		return op; // placeholder
	}

	/**
	 * Performs mean curvature flow on the input mesh with timestep h.
	 * @method module:Projects.MeanCurvatureFlow#integrate
	 * @param {number} h The timestep.
	 */
	integrate(h) {
		let vertices = this.geometry.mesh.vertices;
		
		let M = this.geometry.massMatrix(this.vertexIndex)
		let A = this.buildFlowOperator(M, h)
		let f = DenseMatrix.zeros(vertices.length, 3);
		for(let row = 0; row < vertices.length; row++){
			let pos = this.geometry.positions[row]
			f.set(pos.x, row, 0)
			f.set(pos.y, row, 1)
			f.set(pos.z, row, 2)
		}
		let Mf = M.timesDense(f)
		let out = A.chol().solvePositiveDefinite(Mf)
		//out = M.invertDiagonal().timesDense(out)
		let newVec = {}
		for(let i = 0; i < vertices.length; i++) {
			newVec[i] = new Vector(out.get(i, 0), out.get(i,1), out.get(i,2))
		}
		this.geometry.positions = newVec
		// center mesh positions around origin
		normalize(this.geometry.positions, vertices, false);
	}
}
