"use strict";

class Geometry {
	/**
	 * This class represents the geometry of a {@link module:Core.Mesh Mesh}. This includes information such
	 * as the position of vertices as well as methods to compute edge lengths, corner
	 * angles, face area, normals, discrete curvatures etc.
	 * @constructor module:Core.Geometry
	 * @param {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @param {module:LinearAlgebra.Vector[]} positions An array containing the position of each vertex in a mesh.
	 * @param {boolean} normalizePositions flag to indicate whether positions should be normalized. Default value is true.
	 * @property {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @property {Object} positions A dictionary mapping each vertex to a normalized position.
	 */
	constructor(mesh, positions, normalizePositions = true) {
		this.mesh = mesh;
		this.positions = {};
		for (let i = 0; i < positions.length; i++) {
			let v = this.mesh.vertices[i];
			let p = positions[i];

			this.positions[v] = p;
		}

		if (normalizePositions) {
			normalize(this.positions, mesh.vertices);
		}
	}

	/**
	 * Computes the vector along a halfedge.
	 * @method module:Core.Geometry#vector
	 * @param {module:Core.Halfedge} h The halfedge along which the vector needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vector(h) {
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];

		return b.minus(a);
	}

	/**
	 * Computes the length of an edge.
	 * @method module:Core.Geometry#length
	 * @param {module:Core.Edge} e The edge whose length needs to be computed.
	 * @returns {number}
	 */
	length(e) {
		return this.vector(e.halfedge).norm();
	}

	/**
	 * Computes the midpoint of an edge.
	 * @method module:Core.Geometry#midpoint
	 * @param {module:Core.Edge} e The edge whose midpoint needs to be computed.
	 * @returns {number}
	 */
	midpoint(e) {
		let h = e.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.twin.vertex];

		return (a.plus(b)).over(2);
	}

	/**
	 * Computes the mean edge length of all the edges in a mesh.
	 * @method module:Core.Geometry#meanEdgeLength
	 * @returns {number}
	 */
	meanEdgeLength() {
		let sum = 0;
		let edges = this.mesh.edges;
		for (let e of edges) {
			sum += this.length(e);
		}

		return sum / edges.length;
	}

	/**
	 * Computes the area of a face.
	 * @method module:Core.Geometry#area
	 * @param {module:Core.Face} f The face whose area needs to be computed.
	 * @returns {number}
	 */
	area(f) {
		if (f.isBoundaryLoop()) return 0.0;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return 0.5 * u.cross(v).norm();
	}

	/**
	 * Computes the total surface area of a mesh.
	 * @method module:Core.Geometry#totalArea
	 * @returns {number}
	 */
	totalArea() {
		let sum = 0.0;
		for (let f of this.mesh.faces) {
			sum += this.area(f);
		}

		return sum;
	}

	/**
	 * Computes the normal of a face.
	 * @method module:Core.Geometry#faceNormal
	 * @param {module:Core.Face} f The face whose normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	faceNormal(f) {
		if (f.isBoundaryLoop()) return undefined;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return u.cross(v).unit();
	}

	/**
	 * Computes the centroid of a face.
	 * @method module:Core.Geometry#centroid
	 * @param {module:Core.Face} f The face whose centroid needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	centroid(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		return a.plus(b).plus(c).over(3);
	}

	/**
	 * Computes the circumcenter of a face.
	 * @method module:Core.Geometry#circumcenter
	 * @param {module:Core.Face} f The face whose circumcenter needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	circumcenter(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		let ac = c.minus(a);
		let ab = b.minus(a);
		let w = ab.cross(ac);

		let u = (w.cross(ab)).times(ac.norm2());
		let v = (ac.cross(w)).times(ab.norm2());
		let x = (u.plus(v)).over(2 * w.norm2());

		return x.plus(a);
	}

	/**
	 * Computes an orthonormal bases for a face.
	 * @method module:Core.Geometry#orthonormalBases
	 * @param {module:Core.Face} f The face on which the orthonormal bases needs to be computed.
	 * @returns {module:LinearAlgebra.Vector[]} An array containing two orthonormal vectors tangent to the face.
	 */
	orthonormalBases(f) {
		let e1 = this.vector(f.halfedge).unit();

		let normal = this.faceNormal(f);
		let e2 = normal.cross(e1);

		return [e1, e2];
	}

	/**
	 * Computes the angle (in radians) at a corner.
	 * @method module:Core.Geometry#angle
	 * @param {module:Core.Corner} c The corner at which the angle needs to be computed.
	 * @returns {number} The angle clamped between 0 and π.
	 */
	angle(c) {
		const h = c.halfedge
		const cotan = this.cotan(h)
		const tan = 1.0 / cotan
		let angle = Math.atan(tan)
		if(angle < 0){
			angle = Math.PI + angle
		}
		return angle
	}

	/**
	 * Computes the cotangent of the angle opposite to a halfedge.
	 * @method module:Core.Geometry#cotan
	 * @param {module:Core.Halfedge} h The halfedge opposite to the angle whose cotangent needs to be computed.
	 * @returns {number}
	 */
	cotan(h) {
		if(h.onBoundary){
			return 0
		}
		let A = this.vector(h.next.twin)
		let B = this.vector(h.next.next)

		let cross = A.cross(B)
		let dot = A.dot(B)

		let cotTh = dot / cross.norm() 

		return cotTh // placeholder
	}

	/**
	 * Computes the signed angle (in radians) between two adjacent faces.
	 * @method module:Core.Geometry#dihedralAngle
	 * @param {module:Core.Halfedge} h The halfedge (shared by the two adjacent faces) on which
	 * the dihedral angle is computed.
	 * @returns {number} The dihedral angle.
	 */
	dihedralAngle(h) {
		const a = h.face
		const b = h.twin.face

		const na = this.faceNormal(a)
		const nb = this.faceNormal(b)

		let e = this.vector(h)
		e.normalize()

		const rhs = na.dot(nb)
		const cross = na.cross(nb)
		const lhs = e.dot(cross)
		const angle = Math.atan2(lhs, rhs)
		return angle
	}

	/**
	 * Computes the barycentric dual area of a vertex.
	 * @method module:Core.Geometry#barycentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose barycentric dual area needs to be computed.
	 * @returns {number}
	 */
	barycentricDualArea(v) {
		// TODO
		let faces = [...v.adjacentFaces()]
		let areas = faces.map(f => this.area(f))
		let baryArea = areas.reduce((a,b) => a + b);
		return baryArea/3.0 // placeholder
	}

	/**
	 * Computes the circumcentric dual area of a vertex.
	 * @see {@link http://www.cs.cmu.edu/~kmcrane/Projects/Other/TriangleAreasCheatSheet.pdf}
	 * @method module:Core.Geometry#circumcentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose circumcentric dual area needs to be computed.
	 * @returns {number}
	 */
	circumcentricDualArea(v) {
		let accum = 0.0
		for(let e of v.adjacentEdges()){
			let cotan1 = this.cotan(e.halfedge)
			let cotan2 = this.cotan(e.halfedge.twin)
			let length = this.length(e)
			accum += length * length * cotan1
			accum += length * length * cotan2
		}
		return accum / 8.0; // placeholder
	}

	/**
	 * Computes the normal at a vertex using the "equally weighted" method.
	 * @method module:Core.Geometry#vertexNormalEquallyWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalEquallyWeighted(v) {
		let n = new Vector();
		for (let f of v.adjacentFaces()) {
			let normal = this.faceNormal(f);
			n.incrementBy(normal);
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "face area weights" method.
	 * @method module:Core.Geometry#vertexNormalAreaWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAreaWeighted(v) {
		let n = new Vector()
		for (let f of v.adjacentFaces()){
			const faceArea = this.area(f)
			const normal = this.faceNormal(f)
			normal.scaleBy(faceArea)
			n.incrementBy(normal)
		}
		n.normalize()
		return n
	}

	/**
	 * Computes the normal at a vertex using the "tip angle weights" method.
	 * @method module:Core.Geometry#vertexNormalAngleWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAngleWeighted(v) {
		let n = new Vector();
		let faces = [...v.adjacentFaces()]
		let corners = [...v.adjacentCorners()]
		for (let i = 0; i < faces.length; i++) {
			let f = faces[i]
			let c = corners[i]
			let normal = this.faceNormal(f)
			let angle = this.angle(c)
			normal.scaleBy(angle)
			n.incrementBy(normal)
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "gauss curvature" method.
	 * @method module:Core.Geometry#vertexNormalGaussCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalGaussCurvature(v) {
		let n = new Vector()
		for(let h of v.adjacentHalfedges()) {
			const dihedral = this.dihedralAngle(h)
			let e = this.vector(h)
			e.normalize()
			e.scaleBy(dihedral)
			e.scaleBy(-1)
			n.incrementBy(e)
		}
		n.normalize()
		return n
	}

	/**
	 * Computes the normal at a vertex using the "mean curvature" method (same as the "area gradient" method).
	 * @method module:Core.Geometry#vertexNormalMeanCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalMeanCurvature(v) {
		let n = new Vector()
		for(let h of v.adjacentHalfedges()){
			const cotalpha = this.cotan(h)
			const cotbeta = this.cotan(h.twin)
			let e = this.vector(h)
			e.scaleBy(cotalpha + cotbeta)
			n.incrementBy(e)
		}
		n.normalize()
		return n; // placeholder
	}

	/**
	 * Computes the normal at a vertex using the "inscribed sphere" method.
	 * @method module:Core.Geometry#vertexNormalSphereInscribed
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalSphereInscribed(v) {
		let norm = new Vector()
		let edges = [...v.adjacentHalfedges()]
		for(let i = 0 ; i < edges.length; i++){
			let h = edges[i]
			let nextIdx = (i + 1) % edges.length
			let nexth = edges[nextIdx]

			let hV = this.vector(h)
			let nextHv = this.vector(nexth)
			let cross = nextHv.cross(hV)

			let hVLen2 = hV.norm2()
			let nextHVLen2 = nextHv.norm2()

			cross.scaleBy(1/(hVLen2*nextHVLen2))
			norm.incrementBy(cross)
		}
		norm.normalize()
		return norm
	}

	/**
	 * Computes the angle defect at a vertex (= 2π minus the sum of incident angles
	 * at an interior vertex or π minus the sum of incident angles at a boundary vertex).
	 * @method module:Core.Geometry#angleDefect
	 * @param {module:Core.Vertex} v The vertex whose angle defect needs to be computed.
	 * @returns {number}
	 */
	angleDefect(v) {
		// http://brickisland.net/DDGSpring2019/wp-content/uploads/2019/03/DDG_458_SP19_Lecture15_DiscreteCurvatureI-1.pdf
	
		let angle = 0.0
		let corners = [...v.adjacentCorners()]
		for(let i = 0; i < corners.length; i++){
			const c = corners[i]
			angle += this.angle(c)
		}
		if(v.onBoundary()){
			return Math.PI - angle 
		}
		else{
			return 2*Math.PI - angle
		}
	}

	/**
	 * Computes the (integrated) scalar gauss curvature at a vertex.
	 * @method module:Core.Geometry#scalarGaussCurvature
	 * @param {module:Core.Vertex} v The vertex whose gauss curvature needs to be computed.
	 * @returns {number}
	 */
	scalarGaussCurvature(v) {
		return this.angleDefect(v);
	}

	/**
	 * Computes the (integrated) scalar mean curvature at a vertex.
	 * @method module:Core.Geometry#scalarMeanCurvature
	 * @param {module:Core.Vertex} v The vertex whose mean curvature needs to be computed.
	 * @returns {number}
	 */
	scalarMeanCurvature(v) {
		let total = 0.0
		for(let h of v.adjacentHalfedges()){
			total += this.dihedralAngle(h) * this.length(h.edge)
		}
		return 0.5*total; // placeholder
	}

	/**
	 * Computes the total angle defect (= 2π times the euler characteristic of the mesh).
	 * @method module:Core.Geometry#totalAngleDefect
	 * @returns {number}
	 */
	totalAngleDefect() {

		const total = this.mesh.vertices.map((v) => this.angleDefect(v))
			.reduce((a,b) => a+b)
		return total; // placeholder
	}

	/**
	 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
	 * @method module:Core.Geometry#principalCurvatures
	 * @param {module:Core.Vertex} v The vertex on which the principal curvatures need to be computed.
	 * @returns {number[]} An array containing the minimum and maximum principal curvature values at a vertex.
	 */
	principalCurvatures(v) {
		let H = this.scalarMeanCurvature(v) / this.circumcentricDualArea(v)
		let K = this.scalarGaussCurvature(v) / this.circumcentricDualArea(v)
		let det =  Math.sqrt(H*H - K)
		let sol1 = H + det
		let check1 = K / sol1
		let sol2 = H - det
		let check2 = K / sol2
		let answer = [sol1,sol2]
		let checkH = (sol1+sol2)*.5
		let checkK = sol1*sol2
		answer.sort((a,b) => Math.abs(a) - Math.abs(b))
		return answer // placeholder
	}

	/**
	 * Builds a sparse laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#laplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	laplaceMatrix(vertexIndex) {

		let vertices = this.mesh.vertices
		let T = new Triplet(vertices.length, vertices.length)
		for(let v of vertices){
			let diag = 0
			const row = vertexIndex[v]
			for(let h of v.adjacentHalfedges()){
				let alpha = this.cotan(h)
				let beta = this.cotan(h.twin)
				const col = vertexIndex[h.twin.vertex]
				T.addEntry(-.5*(alpha+beta), row, col)
				diag += .5*(alpha+beta)
			}
			T.addEntry(diag, row, row)
		}
		let mat = SparseMatrix.fromTriplet(T)
		this.perturbMatrixDiagonal(mat)
		return mat
	}

	perturbMatrixDiagonal(matrix) {
		let n = Math.min(matrix.nRows(), matrix.nCols())
		let diagElems = DenseMatrix.ones(n, 1)
		diagElems.scaleBy(1e-8)
		let perturb = SparseMatrix.diag(diagElems)
		matrix.incrementBy(perturb)
	}

	/**
	 * Builds a sparse diagonal mass matrix containing the barycentric dual area of each vertex
	 * of a mesh.
	 * @method module:Core.Geometry#massMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	massMatrix(vertexIndex) {
		let vertices = this.mesh.vertices
		let T = new Triplet(vertices.length, vertices.length)

		vertices.forEach(v => {
			let bary = this.barycentricDualArea(v)
			T.addEntry(bary, vertexIndex[v], vertexIndex[v])
		});

		let adjMat = SparseMatrix.fromTriplet(T)
		return adjMat
	}

	/**
	 * Builds a sparse complex laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#complexLaplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	complexLaplaceMatrix(vertexIndex) {
		let vertices = this.mesh.vertices
		let T = new ComplexTriplet(vertices.length, vertices.length)
		for(let v of vertices){
			let diag = new Complex(0, 0)
			const row = vertexIndex[v]
			for(let h of v.adjacentHalfedges()){
				let alpha = this.cotan(h)
				let beta = this.cotan(h.twin)
				const col = vertexIndex[h.twin.vertex]
				let insert = new Complex(.25*(alpha+beta), 0)
				let negInsert = insert.timesReal(-1)
				T.addEntry(negInsert, row, col)
				diag = diag.plus(insert)
			}
			T.addEntry(diag, row, row)
		}
		let mat = ComplexSparseMatrix.fromTriplet(T)
		this.complexPerturbMatrixDiagonal(mat)
		return mat
	}

	complexPerturbMatrixDiagonal(matrix) {
		let n = Math.min(matrix.nRows(), matrix.nCols())
		let diagElems = ComplexDenseMatrix.ones(n, 1)
		diagElems.scaleBy(new Complex(1e-8, 0))
		let perturb = ComplexSparseMatrix.diag(diagElems)
		matrix.incrementBy(perturb)
	}
}

/**
 * Centers a mesh about the origin and rescales it to unit radius.
 * @global
 * @function module:Core.normalize
 * @param {module:LinearAlgebra.Vector[]} positions The position of each vertex in the vertices array.
 * @param {module:Core.Vertex[]} vertices The vertices of a mesh.
 * @param {boolean} rescale A flag indicating whether mesh positions should be scaled to a unit radius.
 */
function normalize(positions, vertices, rescale = true) {
	// compute center of mass
	let N = vertices.length;
	let cm = new Vector();
	for (let v of vertices) {
		let p = positions[v];

		cm.incrementBy(p);
	}
	cm.divideBy(N);

	// translate to origin and determine radius
	let radius = -1;
	for (let v of vertices) {
		let p = positions[v];

		p.decrementBy(cm);
		radius = Math.max(radius, p.norm());
	}

	// rescale to unit radius
	if (rescale) {
		for (let v of vertices) {
			let p = positions[v];

			p.divideBy(radius);
		}
	}
}