"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		// TODO
		let massMat = geometry.massMatrix()
		return massMat
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		
		const edges = geometry.mesh.edges;
		let T = new Triplet(edges.length, edges.length)
		edges.forEach((edge) => {
			const h = edge.halfedge;
			const alpha = geometry.cotan(h)
			const beta = geometry.cotan(h.twin)
			const ratio = (alpha + beta) / 2
			T.addEntry(ratio, edgeIndex[edge], edgeIndex[edge])
		})
		const dual = SparseMatrix.fromTriplet(T)
		return dual
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		const faces = geometry.mesh.faces
		let T = new Triplet(faces.length, faces.length)
		
		for(let f of faces){
			let geomArea = geometry.area(f)
			T.addEntry(1/geomArea, faceIndex[f], faceIndex[f])
		}
		return SparseMatrix.fromTriplet(T)
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		// TODO
		const edges = geometry.mesh.edges
		const vertices = geometry.mesh.vertices
		let T = new Triplet(edges.length, vertices.length)
		for(let e of edges){
			let v1 = e.halfedge.vertex
			let v2 = e.halfedge.twin.vertex
			let row = edgeIndex[e]
			// alpha(v1) - alpha(v2)
			T.addEntry( 1, row, vertexIndex[v1])
			T.addEntry(-1, row, vertexIndex[v2])
		}
		return SparseMatrix.fromTriplet(T); // placeholder
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		const edges = geometry.mesh.edges
		const faces = geometry.mesh.faces
		let T = new Triplet(faces.length, edges.length)
		for(let f of faces){
			const row = faceIndex[f]
			for(let e of f.adjacentEdges()){
				// get edge orientation here. Not sure if I understand this fully.
				// see http://brickisland.net/DDGSpring2019/2019/02/13/1669/#comment-107
				const entry = e.halfedge.face == f? 1 : -1
				T.addEntry(entry, row, edgeIndex[e])
			}
		}
		return SparseMatrix.fromTriplet(T)
	}
}
