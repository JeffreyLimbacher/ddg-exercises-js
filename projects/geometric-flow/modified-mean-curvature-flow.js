"use strict";

class ModifiedMeanCurvatureFlow extends MeanCurvatureFlow {
	/**
	 * This class performs a {@link http://cs.jhu.edu/~misha/MyPapers/SGP12.pdf modified version} of {@link https://www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf mean curvature flow} on a surface mesh.
	 * @constructor module:Projects.ModifiedMeanCurvatureFlow
	 * @augments module:Projects.MeanCurvatureFlow
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 */
	constructor(geometry) {
		super(geometry);

		this.vertexIndex = indexElements(geometry.mesh.vertices);
		this.A = geometry.laplaceMatrix(this.vertexIndex);
	}

	/**
	 * @inheritdoc
	 */
	buildFlowOperator(M, h) {
		// TODO
		let L = this.A.timesReal(h)
		let op = M.plus(L)
		return op; // placeholder
	}
}
