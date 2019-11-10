"use strict";

class SpectralConformalParameterization {
	/**
	 * This class implements the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf spectral conformal parameterization} algorithm to flatten
	 * surface meshes with boundaries conformally.
	 * @constructor module:Projects.SpectralConformalParameterization
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the complex conformal energy matrix EC = ED - A.
	 * @private
	 * @method module:Projects.SpectralConformalParameterization#buildConformalEnergy
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	buildConformalEnergy() {
		let L = this.geometry.complexLaplaceMatrix(this.vertexIndex)
		let n = this.geometry.mesh.vertices.length
		let A = new ComplexTriplet(n, n)
		for(let bFace of this.geometry.mesh.boundaries){
			for(let e of bFace.adjacentEdges()){
				let h1 = e.halfedge
				let h2 = h1.twin
				let v1 = h1.vertex
				let v2 = h2.vertex
				let vi1 = this.vertexIndex[v1]
				let vi2 = this.vertexIndex[v2]
				A.addEntry(new Complex(0, -1/4), vi1, vi2)
				A.addEntry(new Complex(0, 1/4), vi2, vi1)
			}
		}
		let Ac = ComplexSparseMatrix.fromTriplet(A)
		L.decrementBy(Ac)
		return L
	}

	/**
	 * Flattens the input surface mesh with 1 or more boundaries conformally.
	 * @method module:Projects.SpectralConformalParameterization#flatten
	 * @returns {Object} A dictionary mapping each vertex to a vector of planar coordinates.
	 */
	flatten() {
		let vertices = this.geometry.mesh.vertices;
		let flattening = this.geometry.positions 
		let x = Solvers.solveInversePowerMethod(this.buildConformalEnergy())
		flattening = {}
		for(let v of this.geometry.mesh.vertices){
			let elem = x.get(this.vertexIndex[v], 0)

			flattening[v] = new Vector(elem.re, elem.im)
		}
		// normalize flattening
		normalize(flattening, vertices);

		return flattening;
	}
}
