"use strict";

/**
 * @module Projects
 */

function eqSet(as, bs) {
        // https://stackoverflow.com/questions/31128855/comparing-ecma6-sets-for-equality
        if (as.size !== bs.size) return false;
        for (var a of as) if (!bs.has(a)) return false;
        return true;
}

class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                //Solution to problem 1
                let arraysToAssign = [mesh.vertices, mesh.edges, mesh.faces];
                arraysToAssign.forEach( (arr) => {
                for(let i = 0; i < arr.length; i++){
                                arr[i].index = i;
                        }
                });
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.DenseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                //Solution to problem 2
                let T = new Triplet(mesh.edges.length, mesh.vertices.length)
                mesh.edges.forEach((edge) => {
                        let row = edge.index;
                        let col1 = edge.halfedge.vertex.index
                        let col2 = edge.halfedge.twin.vertex.index
                        T.addEntry(1, row, col1)
                        T.addEntry(1, row, col2)
                });
                let adjMat = SparseMatrix.fromTriplet(T)
                return adjMat
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.DenseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                let T = new Triplet(mesh.faces.length, mesh.edges.length)
                for(let f of mesh.faces){
                        let row = f.index
                        for(let e of f.adjacentEdges()){
                                let col = e.index
                                T.addEntry(1.0, row, col)
                        }
                }
                let adjMat = SparseMatrix.fromTriplet(T)
                return adjMat
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                let vecSize = this.totalVertices()
                return this.vectorFromSet(subset.vertices, vecSize)
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                let vecSize = this.totalEdges()
                return this.vectorFromSet(subset.edges, vecSize)
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                let vecSize = this.totalFaces()
                return this.vectorFromSet(subset.faces, vecSize)
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                // TODO
                let vertexVec = this.buildVertexVector(subset)
                let edgeVec = this.buildEdgeVector(subset)
                let faceVec = this.buildFaceVector(subset)

                let edgesFromVerts = this.A0.timesSparse(vertexVec)
                let facesFromVerts = this.A1.timesSparse(edgesFromVerts)

                let facesFromEdges = this.A1.timesSparse(edgeVec)

                let allEdges = edgesFromVerts.plus(edgeVec)
                let allFaces = facesFromEdges.plus(facesFromVerts).plus(faceVec)


                let vertsD = vertexVec.toDense()
                let edgesD = allEdges.toDense()
                let facesD = allFaces.toDense()
                
                return this.meshSubsetFromVectors(vertsD, edgesD, facesD)
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                // TODO
                let vertexVec = this.buildVertexVector(subset)
                let edgeVec = this.buildEdgeVector(subset)
                let faceVec = this.buildFaceVector(subset)

                let edgesFromFaces = this.A1.transpose().timesSparse(faceVec)

                let allEdges = edgeVec.plus(edgesFromFaces) 
                let verticesFromEdges = this.A0.transpose().timesSparse(allEdges)
                let allVertices = verticesFromEdges.plus(vertexVec)

                let vertsD = allVertices.toDense()
                let edgesD = allEdges.toDense()
                let facesD = faceVec.toDense()

                return this.meshSubsetFromVectors(vertsD, edgesD, facesD)
        }


        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                // TODO
                let closureSubset = this.closure(subset);
                let starSubset = this.star(subset);

                let closureStar = this.closure(starSubset);
                let starClosure = this.star(closureSubset);

                closureStar.deleteSubset(starClosure);

                return closureStar; // placeholder
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                // The closure of a simplical complex should be itself
                let smallestComplex = this.closure(subset)
                return smallestComplex.equals(subset)
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {

                // A pure complex must be a complex
                if(!this.isComplex(subset)){
                        return -1
                }

                let degree = this.getMaxDegree(subset)
                if( degree >= 2){
                        // Check to make sure every face vector has its corresponding edges
                        let faceVec = this.buildFaceVector(subset)
                        let expectedEdges = this.A1.transpose().timesSparse(faceVec)
                        let expectedEdgesSet = this.setFromVector(expectedEdges.toDense())
                        if(!eqSet(expectedEdgesSet, subset.edges)){
                                return -1
                        }
                }

                if (degree >= 1) {
                        // Check to see if every edge has its corresponding vertices.
                        let edgeVec = this.buildEdgeVector(subset)
                        let expectedVertices = this.A0.transpose().timesSparse(edgeVec)
                        let expectedVerticesSet = this.setFromVector(expectedVertices.toDense())
                        if(!eqSet(expectedVerticesSet, subset.vertices)){
                                return -1
                        }
                }
                return degree
        }

        getMaxDegree(subset){
                if(subset.faces.size > 0){
                        return 2
                } 
                else if (subset.edges.size > 0) {
                        return 1
                }
                else{
                        return 0
                }
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                let degree = this.isPureComplex(subset)
                if(degree === -1){
                        return subset
                }
                // boundary, bd(K), is the closure of the set of all simplices sigma tht are proper faces of exactly one simplex of K'
                // depends on degree (I think) 

                // Technically we can have the same code for each degree (e.g., in a pure 2-complex, no vertex appears in exactly one edge
                // so  the second if branch will return nothing) but we can avoid unnecessary operations and capture the essence of the
                // problem like this, so whatever.

                if(degree === 2){
                        let faceVec = this.buildFaceVector(subset)
                        let edges = this.A1.transpose().timesSparse(faceVec).toDense()
                        // Any element with one as an entry was in exactly one face
                        let boundaryEdges = this.vectorToSetMatchingVal(edges, 1.0)
                        return this.closure(new MeshSubset(new Set(), boundaryEdges, new Set()))
                }
                else if (degree === 1){
                        let edgeVec = this.buildEdgeVector(subset)
                        let vertices = this.A0.transpose().timesSparse(edgeVec).toDense()
                        let boundaryVertices = this.vectorToSetMatchingVal(vertices, 1.0)
                        return new MeshSubset(boundaryVertices, new Set(), new Set())
                }

                // closure of empty set (proper faces of vertices) is empty set
                return new MeshSubset(); 
        }

        vectorToSetMatchingVal(vector, val) {
                let valSet = new Set()
                for(let i = 0; i < vector.nRows(); i++){
                        let vecVal = vector.get(i, 0)
                        if(vecVal == val) valSet.add(i)
                }
                return valSet
        }

        /**
         * Returns the total number of vertices in the mesh.
         * @method module:Projects.SimplicialComplexOperators#totalVertices
         * @returns {Number} The total number of vertices.
         */
        totalVertices() {
                return this.A0.nCols();
        }

        /**
         * Returns the total number of edges in the mesh.
         * @method module:Projects.SimplicialComplexOperators#totalEdges
         * @returns {Number} The total number of edges.
         */
        totalEdges() {
                return this.A0.nRows();
        }

        /**
         * Returns the total number of faces in the mesh.
         * @method module:Projects.SimplicialComplexOperators#totalFaces
         * @returns {Number} The total number of faces.
         */
        totalFaces() {
                return this.A1.nRows();
        }

        meshSubsetFromVectors(vertsD, edgesD, facesD) {
                let verticesSet = this.setFromVector(vertsD)
                let edgesSet = this.setFromVector(edgesD)
                let facesSet = this.setFromVector(facesD)
                return new MeshSubset(verticesSet, edgesSet, facesSet);
        }

        setFromVector(vector){
                //assume column vector
                let range = vector.nRows();
                let vecSet = new Set()
                for(let i = 0; i < range; i++){
                        let val = vector.get(i,0)
                        if(val > 1e-5) vecSet.add(i)
                }
                return vecSet
        }

        vectorFromSet(set, vecSize){
                let T = new Triplet(vecSize, 1)
                for(let v of set){
                        T.addEntry(1.0, v, 0)
                }
                let vertVector = SparseMatrix.fromTriplet(T)
                return vertVector
        }

}
