#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/RadialBasisFctMapping.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;
using namespace precice::mapping;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(RadialBasisFunctionMapping)

void addGlobalIndex(mesh::PtrMesh &mesh, int offset = 0)
{
  for (mesh::Vertex &v : mesh->vertices()) {
    v.setGlobalIndex(v.getID() + offset);
  }
}

BOOST_AUTO_TEST_SUITE(Parallel)

/// Holds rank, owner, position and value of a single vertex
struct VertexSpecification {
  int                 rank;
  int                 owner;
  std::vector<double> position;
  std::vector<double> value;
};

/*
MeshSpecification format:
{ {rank, owner rank, {x, y, z}, {v}}, ... }

also see struct VertexSpecification.

- -1 on rank means all ranks
- -1 on owner rank means no rank
- x, y, z is position of vertex, z is optional, 2D mesh will be created then
- v is the value of the respective vertex. Only 1D supported at this time.

ReferenceSpecification format:
{ {rank, {v}, ... }
- -1 on rank means all ranks
- v is the expected value of n-th vertex on that particular rank
*/
using MeshSpecification = std::vector<VertexSpecification>;

/// Contains which values are expected on which rank: rank -> vector of data.
using ReferenceSpecification = std::vector<std::pair<int, std::vector<double>>>;

void getDistributedMesh(const TestContext &      context,
                        MeshSpecification const &vertices,
                        mesh::PtrMesh &          mesh,
                        mesh::PtrData &          data,
                        int                      globalIndexOffset = 0)
{
  Eigen::VectorXd d;

  int i = 0;
  for (auto &vertex : vertices) {
    if (vertex.rank == context.rank or vertex.rank == -1) {
      if (vertex.position.size() == 3) // 3-dimensional
        mesh->createVertex(Eigen::Vector3d(vertex.position.data()));
      else if (vertex.position.size() == 2) // 2-dimensional
        mesh->createVertex(Eigen::Vector2d(vertex.position.data()));

      int valueDimension = vertex.value.size();

      if (vertex.owner == context.rank)
        mesh->vertices().back().setOwner(true);
      else
        mesh->vertices().back().setOwner(false);

      d.conservativeResize(i * valueDimension + valueDimension);
      // Get data in every dimension
      for (int dim = 0; dim < valueDimension; ++dim) {
        d(i * valueDimension + dim) = vertex.value.at(dim);
      }
      i++;
    }
  }
  addGlobalIndex(mesh, globalIndexOffset);
  mesh->allocateDataValues();
  data->values() = d;
}

void testDistributed(const TestContext &    context,
                     Mapping &              mapping,
                     MeshSpecification      inMeshSpec,
                     MeshSpecification      outMeshSpec,
                     ReferenceSpecification referenceSpec,
                     int                    inGlobalIndexOffset = 0)
{
  int meshDimension  = inMeshSpec.at(0).position.size();
  int valueDimension = inMeshSpec.at(0).value.size();

  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", meshDimension, false, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", valueDimension);
  int           inDataID = inData->getID();

  getDistributedMesh(context, inMeshSpec, inMesh, inData, inGlobalIndexOffset);

  mesh::PtrMesh outMesh(new mesh::Mesh("outMesh", meshDimension, false, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", valueDimension);
  int           outDataID = outData->getID();

  getDistributedMesh(context, outMeshSpec, outMesh, outData);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  mapping.map(inDataID, outDataID);

  int index = 0;
  for (auto &referenceVertex : referenceSpec) {
    if (referenceVertex.first == context.rank or referenceVertex.first == -1) {
      for (int dim = 0; dim < valueDimension; ++dim) {
        BOOST_TEST_INFO("Index of vertex: " << index << " - Dimension: " << dim);
        BOOST_TEST(outData->values()(index * valueDimension + dim) == referenceVertex.second.at(dim));
      }
      ++index;
    }
  }
  BOOST_TEST(outData->values().size() == index * valueDimension);
}

/// Test with a homogenous distribution of mesh amoung ranks
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV1)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);

  testDistributed(context, mapping,
                  {// Consistent mapping: The inMesh is communicated
                   {-1, 0, {0, 0}, {1}},
                   {-1, 0, {0, 1}, {2}},
                   {-1, 1, {1, 0}, {3}},
                   {-1, 1, {1, 1}, {4}},
                   {-1, 2, {2, 0}, {5}},
                   {-1, 2, {2, 1}, {6}},
                   {-1, 3, {3, 0}, {7}},
                   {-1, 3, {3, 1}, {8}}},
                  {// The outMesh is local, distributed amoung all ranks
                   {0, -1, {0, 0}, {0}},
                   {0, -1, {0, 1}, {0}},
                   {1, -1, {1, 0}, {0}},
                   {1, -1, {1, 1}, {0}},
                   {2, -1, {2, 0}, {0}},
                   {2, -1, {2, 1}, {0}},
                   {3, -1, {3, 0}, {0}},
                   {3, -1, {3, 1}, {0}}},
                  {// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                   {0, {1}},
                   {0, {2}},
                   {1, {3}},
                   {1, {4}},
                   {2, {5}},
                   {2, {6}},
                   {3, {7}},
                   {3, {8}}});
}

BOOST_AUTO_TEST_CASE(DistributedConsistent2DV1Vector)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);

  testDistributed(context, mapping,
                  {// Consistent mapping: The inMesh is communicated
                   {-1, 0, {0, 0}, {1, 4}},
                   {-1, 0, {0, 1}, {2, 5}},
                   {-1, 1, {1, 0}, {3, 6}},
                   {-1, 1, {1, 1}, {4, 7}},
                   {-1, 2, {2, 0}, {5, 8}},
                   {-1, 2, {2, 1}, {6, 9}},
                   {-1, 3, {3, 0}, {7, 10}},
                   {-1, 3, {3, 1}, {8, 11}}},
                  {// The outMesh is local, distributed amoung all ranks
                   {0, -1, {0, 0}, {0, 0}},
                   {0, -1, {0, 1}, {0, 0}},
                   {1, -1, {1, 0}, {0, 0}},
                   {1, -1, {1, 1}, {0, 0}},
                   {2, -1, {2, 0}, {0, 0}},
                   {2, -1, {2, 1}, {0, 0}},
                   {3, -1, {3, 0}, {0, 0}},
                   {3, -1, {3, 1}, {0, 0}}},
                  {// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                   {0, {1, 4}},
                   {0, {2, 5}},
                   {1, {3, 6}},
                   {1, {4, 7}},
                   {2, {5, 8}},
                   {2, {6, 9}},
                   {3, {7, 10}},
                   {3, {8, 11}}});
}

/// Using a more heterogenous distributon of vertices and owner
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV2)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);

  testDistributed(context, mapping,
                  {// Consistent mapping: The inMesh is communicated, rank 2 owns no vertices
                   {-1, 0, {0, 0}, {1}},
                   {-1, 0, {0, 1}, {2}},
                   {-1, 1, {1, 0}, {3}},
                   {-1, 1, {1, 1}, {4}},
                   {-1, 1, {2, 0}, {5}},
                   {-1, 3, {2, 1}, {6}},
                   {-1, 3, {3, 0}, {7}},
                   {-1, 3, {3, 1}, {8}}},
                  {// The outMesh is local, rank 1 is empty
                   {0, -1, {0, 0}, {0}},
                   {0, -1, {0, 1}, {0}},
                   {0, -1, {1, 0}, {0}},
                   {2, -1, {1, 1}, {0}},
                   {2, -1, {2, 0}, {0}},
                   {2, -1, {2, 1}, {0}},
                   {3, -1, {3, 0}, {0}},
                   {3, -1, {3, 1}, {0}}},
                  {// Tests for {0, 1, 2} on the first rank,
                   // second rank (consistent with the outMesh) is empty, ...
                   {0, {1}},
                   {0, {2}},
                   {0, {3}},
                   {2, {4}},
                   {2, {5}},
                   {2, {6}},
                   {3, {7}},
                   {3, {8}}});
}

/// Test with a very heterogenous distributed and non-continues ownership
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV3)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 0, 0, 4};

  testDistributed(context, mapping,
                  {
                      // Rank 0 has part of the mesh, owns a subpart
                      {0, 0, {0, 0}, {1}},
                      {0, 0, {0, 1}, {2}},
                      {0, 0, {1, 0}, {3}},
                      {0, -1, {1, 1}, {4}},
                      {0, -1, {2, 0}, {5}},
                      {0, -1, {2, 1}, {6}},
                      // Rank 1 has no vertices
                      // Rank 2 has the entire mesh, but owns just 3 and 5.
                      {2, -1, {0, 0}, {1}},
                      {2, -1, {0, 1}, {2}},
                      {2, -1, {1, 0}, {3}},
                      {2, 2, {1, 1}, {4}},
                      {2, -1, {2, 0}, {5}},
                      {2, 2, {2, 1}, {6}},
                      {2, -1, {3, 0}, {7}},
                      {2, -1, {3, 1}, {8}},
                      // Rank 3 has the last 4 vertices, owns 4, 6 and 7
                      {3, 3, {2, 0}, {5}},
                      {3, -1, {2, 1}, {6}},
                      {3, 3, {3, 0}, {7}},
                      {3, 3, {3, 1}, {8}},
                  },
                  {// The outMesh is local, rank 1 is empty
                   {0, -1, {0, 0}, {0}},
                   {0, -1, {0, 1}, {0}},
                   {0, -1, {1, 0}, {0}},
                   {2, -1, {1, 1}, {0}},
                   {2, -1, {2, 0}, {0}},
                   {2, -1, {2, 1}, {0}},
                   {3, -1, {3, 0}, {0}},
                   {3, -1, {3, 1}, {0}}},
                  {// Tests for {0, 1, 2} on the first rank,
                   // second rank (consistent with the outMesh) is empty, ...
                   {0, {1}},
                   {0, {2}},
                   {0, {3}},
                   {2, {4}},
                   {2, {5}},
                   {2, {6}},
                   {3, {7}},
                   {3, {8}}},
                  globalIndexOffsets.at(context.rank));
}

/// Test with a very heterogenous distributed and non-continues ownership
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV3Vector)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 0, 0, 4};

  testDistributed(context, mapping,
                  {
                      // Rank 0 has part of the mesh, owns a subpart
                      {0, 0, {0, 0}, {1, 4}},
                      {0, 0, {0, 1}, {2, 5}},
                      {0, 0, {1, 0}, {3, 6}},
                      {0, -1, {1, 1}, {4, 7}},
                      {0, -1, {2, 0}, {5, 8}},
                      {0, -1, {2, 1}, {6, 9}},
                      // Rank 1 has no vertices
                      // Rank 2 has the entire mesh, but owns just 3 and 5.
                      {2, -1, {0, 0}, {1, 4}},
                      {2, -1, {0, 1}, {2, 5}},
                      {2, -1, {1, 0}, {3, 6}},
                      {2, 2, {1, 1}, {4, 7}},
                      {2, -1, {2, 0}, {5, 8}},
                      {2, 2, {2, 1}, {6, 9}},
                      {2, -1, {3, 0}, {7, 10}},
                      {2, -1, {3, 1}, {8, 11}},
                      // Rank 3 has the last 4 vertices, owns 4, 6 and 7
                      {3, 3, {2, 0}, {5, 8}},
                      {3, -1, {2, 1}, {6, 9}},
                      {3, 3, {3, 0}, {7, 10}},
                      {3, 3, {3, 1}, {8, 11}},
                  },
                  {// The outMesh is local, rank 1 is empty
                   {0, -1, {0, 0}, {0, 0}},
                   {0, -1, {0, 1}, {0, 0}},
                   {0, -1, {1, 0}, {0, 0}},
                   {2, -1, {1, 1}, {0, 0}},
                   {2, -1, {2, 0}, {0, 0}},
                   {2, -1, {2, 1}, {0, 0}},
                   {3, -1, {3, 0}, {0, 0}},
                   {3, -1, {3, 1}, {0, 0}}},
                  {// Tests for {0, 1, 2} on the first rank,
                   // second rank (consistent with the outMesh) is empty, ...
                   {0, {1, 4}},
                   {0, {2, 5}},
                   {0, {3, 6}},
                   {2, {4, 7}},
                   {2, {5, 8}},
                   {2, {6, 9}},
                   {3, {7, 10}},
                   {3, {8, 11}}},
                  globalIndexOffsets.at(context.rank));
}

/// Some ranks are empty, does not converge
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV4)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  ThinPlateSplines                        fct;
  RadialBasisFctMapping<ThinPlateSplines> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 0, 0, 0};

  testDistributed(context, mapping,
                  {
                      // Rank 0 has no vertices
                      // Rank 1 has the entire mesh, owns a subpart
                      {1, 1, {0, 0}, {1.1}},
                      {1, 1, {0, 1}, {2.5}},
                      {1, 1, {1, 0}, {3}},
                      {1, 1, {1, 1}, {4}},
                      {1, -1, {2, 0}, {5}},
                      {1, -1, {2, 1}, {6}},
                      {1, -1, {3, 0}, {7}},
                      {1, -1, {3, 1}, {8}},
                      // Rank 2 has the entire mesh, owns a subpart
                      {2, -1, {0, 0}, {1.1}},
                      {2, -1, {0, 1}, {2.5}},
                      {2, -1, {1, 0}, {3}},
                      {2, -1, {1, 1}, {4}},
                      {2, 2, {2, 0}, {5}},
                      {2, 2, {2, 1}, {6}},
                      {2, 2, {3, 0}, {7}},
                      {2, 2, {3, 1}, {8}},
                      // Rank 3 has no vertices
                  },
                  {// The outMesh is local, rank 0 and 3 are empty
                   // not same order as input mesh and vertex (2,0) appears twice
                   {1, -1, {2, 0}, {0}},
                   {1, -1, {1, 0}, {0}},
                   {1, -1, {0, 1}, {0}},
                   {1, -1, {1, 1}, {0}},
                   {1, -1, {0, 0}, {0}},
                   {2, -1, {2, 0}, {0}},
                   {2, -1, {2, 1}, {0}},
                   {2, -1, {3, 0}, {0}},
                   {2, -1, {3, 1}, {0}}},
                  {{1, {5}},
                   {1, {3}},
                   {1, {2.5}},
                   {1, {4}},
                   {1, {1.1}},
                   {2, {5}},
                   {2, {6}},
                   {2, {7}},
                   {2, {8}}},
                  globalIndexOffsets.at(context.rank));
}

// same as 2DV4, but all ranks have vertices
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV5)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  ThinPlateSplines                        fct;
  RadialBasisFctMapping<ThinPlateSplines> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 0, 0, 0};

  testDistributed(context, mapping,
                  {
                      // Every rank has the entire mesh and owns a subpart
                      {0, 0, {0, 0}, {1.1}},
                      {0, 0, {0, 1}, {2.5}},
                      {0, -1, {1, 0}, {3}},
                      {0, -1, {1, 1}, {4}},
                      {0, -1, {2, 0}, {5}},
                      {0, -1, {2, 1}, {6}},
                      {0, -1, {3, 0}, {7}},
                      {0, -1, {3, 1}, {8}},
                      {1, -1, {0, 0}, {1.1}},
                      {1, -1, {0, 1}, {2.5}},
                      {1, 1, {1, 0}, {3}},
                      {1, 1, {1, 1}, {4}},
                      {1, -1, {2, 0}, {5}},
                      {1, -1, {2, 1}, {6}},
                      {1, -1, {3, 0}, {7}},
                      {1, -1, {3, 1}, {8}},
                      {2, -1, {0, 0}, {1.1}},
                      {2, -1, {0, 1}, {2.5}},
                      {2, -1, {1, 0}, {3}},
                      {2, -1, {1, 1}, {4}},
                      {2, 2, {2, 0}, {5}},
                      {2, 2, {2, 1}, {6}},
                      {2, -1, {3, 0}, {7}},
                      {2, -1, {3, 1}, {8}},
                      {3, -1, {0, 0}, {1.1}},
                      {3, -1, {0, 1}, {2.5}},
                      {3, -1, {1, 0}, {3}},
                      {3, -1, {1, 1}, {4}},
                      {3, -1, {2, 0}, {5}},
                      {3, -1, {2, 1}, {6}},
                      {3, 3, {3, 0}, {7}},
                      {3, 3, {3, 1}, {8}},
                  },
                  {// The outMesh is local, rank 0 and 3 are empty
                   // not same order as input mesh and vertex (2,0) appears twice
                   {1, -1, {2, 0}, {0}},
                   {1, -1, {1, 0}, {0}},
                   {1, -1, {0, 1}, {0}},
                   {1, -1, {1, 1}, {0}},
                   {1, -1, {0, 0}, {0}},
                   {2, -1, {2, 0}, {0}},
                   {2, -1, {2, 1}, {0}},
                   {2, -1, {3, 0}, {0}},
                   {2, -1, {3, 1}, {0}}},
                  {{1, {5}},
                   {1, {3}},
                   {1, {2.5}},
                   {1, {4}},
                   {1, {1.1}},
                   {2, {5}},
                   {2, {6}},
                   {2, {7}},
                   {2, {8}}},
                  globalIndexOffsets.at(context.rank));
}

/// same as 2DV4, but strictly linear input values, converges and gives correct results
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV6,
                     *boost::unit_test::tolerance(1e-7))
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  ThinPlateSplines                        fct;
  RadialBasisFctMapping<ThinPlateSplines> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 0, 0, 0};

  testDistributed(context, mapping,
                  {
                      // Rank 0 has no vertices
                      // Rank 1 has the entire mesh, owns a subpart
                      {1, 1, {0, 0}, {1}},
                      {1, 1, {0, 1}, {2}},
                      {1, 1, {1, 0}, {3}},
                      {1, 1, {1, 1}, {4}},
                      {1, -1, {2, 0}, {5}},
                      {1, -1, {2, 1}, {6}},
                      {1, -1, {3, 0}, {7}},
                      {1, -1, {3, 1}, {8}},
                      // Rank 2 has the entire mesh, owns a subpart
                      {2, -1, {0, 0}, {1}},
                      {2, -1, {0, 1}, {2}},
                      {2, -1, {1, 0}, {3}},
                      {2, -1, {1, 1}, {4}},
                      {2, 2, {2, 0}, {5}},
                      {2, 2, {2, 1}, {6}},
                      {2, 2, {3, 0}, {7}},
                      {2, 2, {3, 1}, {8}},
                      // Rank 3 has no vertices
                  },
                  {// The outMesh is local, rank 0 and 3 are empty
                   // not same order as input mesh and vertex (2,0) appears twice
                   {1, -1, {2, 0}, {0}},
                   {1, -1, {1, 0}, {0}},
                   {1, -1, {0, 1}, {0}},
                   {1, -1, {1, 1}, {0}},
                   {1, -1, {0, 0}, {0}},
                   {2, -1, {2, 0}, {0}},
                   {2, -1, {2, 1}, {0}},
                   {2, -1, {3, 0}, {0}},
                   {2, -1, {3, 1}, {0}}},
                  {{1, {5}},
                   {1, {3}},
                   {1, {2}},
                   {1, {4}},
                   {1, {1}},
                   {2, {5}},
                   {2, {6}},
                   {2, {7}},
                   {2, {8}}},
                  globalIndexOffsets.at(context.rank));
}

/// Test with a homogenous distribution of mesh amoung ranks
BOOST_AUTO_TEST_CASE(DistributedConservative2DV1)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  testDistributed(context, mapping,
                  {// Conservative mapping: The inMesh is local
                   {0, -1, {0, 0}, {1}},
                   {0, -1, {0, 1}, {2}},
                   {1, -1, {1, 0}, {3}},
                   {1, -1, {1, 1}, {4}},
                   {2, -1, {2, 0}, {5}},
                   {2, -1, {2, 1}, {6}},
                   {3, -1, {3, 0}, {7}},
                   {3, -1, {3, 1}, {8}}},
                  {// The outMesh is distributed
                   {-1, 0, {0, 0}, {0}},
                   {-1, 0, {0, 1}, {0}},
                   {-1, 1, {1, 0}, {0}},
                   {-1, 1, {1, 1}, {0}},
                   {-1, 2, {2, 0}, {0}},
                   {-1, 2, {2, 1}, {0}},
                   {-1, 3, {3, 0}, {0}},
                   {-1, 3, {3, 1}, {0}}},
                  {// Tests for {0, 1, 0, 0, 0, 0, 0, 0} on the first rank,
                   // {0, 0, 2, 3, 0, 0, 0, 0} on the second, ...
                   {0, {1}},
                   {0, {2}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {1, {0}},
                   {1, {0}},
                   {1, {3}},
                   {1, {4}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {5}},
                   {2, {6}},
                   {2, {0}},
                   {2, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {7}},
                   {3, {8}}},
                  context.rank * 2);
}

/// Test with a homogenous distribution of mesh amoung ranks
BOOST_AUTO_TEST_CASE(DistributedConservative2DV1Vector)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  testDistributed(context, mapping,
                  {// Conservative mapping: The inMesh is local
                   {0, -1, {0, 0}, {1, 4}},
                   {0, -1, {0, 1}, {2, 5}},
                   {1, -1, {1, 0}, {3, 6}},
                   {1, -1, {1, 1}, {4, 7}},
                   {2, -1, {2, 0}, {5, 8}},
                   {2, -1, {2, 1}, {6, 9}},
                   {3, -1, {3, 0}, {7, 10}},
                   {3, -1, {3, 1}, {8, 11}}},
                  {// The outMesh is distributed
                   {-1, 0, {0, 0}, {0, 0}},
                   {-1, 0, {0, 1}, {0, 0}},
                   {-1, 1, {1, 0}, {0, 0}},
                   {-1, 1, {1, 1}, {0, 0}},
                   {-1, 2, {2, 0}, {0, 0}},
                   {-1, 2, {2, 1}, {0, 0}},
                   {-1, 3, {3, 0}, {0, 0}},
                   {-1, 3, {3, 1}, {0, 0}}},
                  {// Tests for {0, 1, 0, 0, 0, 0, 0, 0} on the first rank,
                   // {0, 0, 2, 3, 0, 0, 0, 0} on the second, ...
                   {0, {1, 4}},
                   {0, {2, 5}},
                   {0, {0, 0}},
                   {0, {0, 0}},
                   {0, {0, 0}},
                   {0, {0, 0}},
                   {0, {0, 0}},
                   {0, {0, 0}},
                   {1, {0, 0}},
                   {1, {0, 0}},
                   {1, {3, 6}},
                   {1, {4, 7}},
                   {1, {0, 0}},
                   {1, {0, 0}},
                   {1, {0, 0}},
                   {1, {0, 0}},
                   {2, {0, 0}},
                   {2, {0, 0}},
                   {2, {0, 0}},
                   {2, {0, 0}},
                   {2, {5, 8}},
                   {2, {6, 9}},
                   {2, {0, 0}},
                   {2, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {7, 10}},
                   {3, {8, 11}}},
                  context.rank * 2);
}

/// Using a more heterogenous distribution of vertices and owner
BOOST_AUTO_TEST_CASE(DistributedConservative2DV2)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves())
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 0, 4, 6};

  testDistributed(context, mapping,
                  {// Conservative mapping: The inMesh is local but rank 0 has no vertices
                   {1, -1, {0, 0}, {1}},
                   {1, -1, {0, 1}, {2}},
                   {1, -1, {1, 0}, {3}},
                   {1, -1, {1, 1}, {4}},
                   {2, -1, {2, 0}, {5}},
                   {2, -1, {2, 1}, {6}},
                   {3, -1, {3, 0}, {7}},
                   {3, -1, {3, 1}, {8}}},
                  {// The outMesh is distributed, rank 0 owns no vertex
                   {-1, 1, {0, 0}, {0}},
                   {-1, 1, {0, 1}, {0}},
                   {-1, 1, {1, 0}, {0}},
                   {-1, 1, {1, 1}, {0}},
                   {-1, 2, {2, 0}, {0}},
                   {-1, 2, {2, 1}, {0}},
                   {-1, 3, {3, 0}, {0}},
                   {-1, 3, {3, 1}, {0}}},
                  {// Tests for {0, 0, 0, 0, 0, 0, 0, 0} on the first rank,
                   // {1, 2, 2, 3, 0, 0, 0, 0} on the second, ...
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {1, {1}},
                   {1, {2}},
                   {1, {3}},
                   {1, {4}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {5}},
                   {2, {6}},
                   {2, {0}},
                   {2, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {7}},
                   {3, {8}}},
                  globalIndexOffsets.at(context.rank));
}

/// Using meshes of different sizes, inMesh is smaller then outMesh
BOOST_AUTO_TEST_CASE(DistributedConservative2DV3)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(2.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 0, 3, 5};

  testDistributed(context, mapping,
                  {// Conservative mapping: The inMesh is local but rank 0 has no vertices
                   {1, -1, {0, 0}, {1}},
                   {1, -1, {1, 0}, {3}},
                   {1, -1, {1, 1}, {4}},
                   {2, -1, {2, 0}, {5}},
                   {2, -1, {2, 1}, {6}},
                   {3, -1, {3, 0}, {7}},
                   {3, -1, {3, 1}, {8}}}, // Sum of all vertices is 34
                  {                       // The outMesh is distributed, rank 0 owns no vertex
                   {-1, 1, {0, 0}, {0}},
                   {-1, 1, {0, 1}, {0}},
                   {-1, 1, {1, 0}, {0}},
                   {-1, 1, {1, 1}, {0}},
                   {-1, 2, {2, 0}, {0}},
                   {-1, 2, {2, 1}, {0}},
                   {-1, 3, {3, 0}, {0}},
                   {-1, 3, {3, 1}, {0}}},
                  {// Tests for {0, 0, 0, 0, 0, 0, 0, 0} on the first rank,
                   // {1, 2, 2, 3, 0, 0, 0, 0} on the second, ...
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {1, {1}},
                   {1, {0}},
                   {1, {3}},
                   {1, {4}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {5}},
                   {2, {6}},
                   {2, {0}},
                   {2, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {7}},
                   {3, {8}}}, // Sum of reference is also 34
                  globalIndexOffsets.at(context.rank));
}

/// Using meshes of different sizes, outMesh is smaller then inMesh
BOOST_AUTO_TEST_CASE(DistributedConservative2DV4,
                     *boost::unit_test::tolerance(1e-6))
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(4.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 2, 4, 6};

  testDistributed(context, mapping,
                  {// Conservative mapping: The inMesh is local
                   {0, -1, {0, 0}, {1}},
                   {0, -1, {0, 1}, {2}},
                   {1, -1, {1, 0}, {3}},
                   {1, -1, {1, 1}, {4}},
                   {2, -1, {2, 0}, {5}},
                   {2, -1, {2, 1}, {6}},
                   {3, -1, {3, 0}, {7}},
                   {3, -1, {3, 1}, {8}}}, // Sum is 36
                  {                       // The outMesh is distributed, rank 0 has no vertex at all
                   {-1, 1, {0, 1}, {0}},
                   {-1, 1, {1, 0}, {0}},
                   {-1, 1, {1, 1}, {0}},
                   {-1, 2, {2, 0}, {0}},
                   {-1, 2, {2, 1}, {0}},
                   {-1, 3, {3, 0}, {0}},
                   {-1, 3, {3, 1}, {0}}},
                  {// Tests for {0, 0, 0, 0, 0, 0, 0, 0} on the first rank,
                   // {2, 3, 4, 3, 0, 0, 0, 0} on the second, ...
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {1, {2.4285714526861519}},
                   {1, {3.61905}},
                   {1, {4.14286}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {5.333333295}},
                   {2, {5.85714}},
                   {2, {0}},
                   {2, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {7.047619}},
                   {3, {7.571428}}}, // Sum is ~36
                  globalIndexOffsets.at(context.rank));
}

/// Tests a non-contigous owner distributed at the outMesh
BOOST_AUTO_TEST_CASE(testDistributedConservative2DV5)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  testDistributed(context, mapping,
                  {// Conservative mapping: The inMesh is local
                   {0, -1, {0, 0}, {1}},
                   {0, -1, {0, 1}, {2}},
                   {1, -1, {1, 0}, {3}},
                   {1, -1, {1, 1}, {4}},
                   {2, -1, {2, 0}, {5}},
                   {2, -1, {2, 1}, {6}},
                   {3, -1, {3, 0}, {7}},
                   {3, -1, {3, 1}, {8}}},
                  {// The outMesh is distributed and non-contigous
                   {-1, 0, {0, 0}, {0}},
                   {-1, 1, {0, 1}, {0}},
                   {-1, 1, {1, 0}, {0}},
                   {-1, 0, {1, 1}, {0}},
                   {-1, 2, {2, 0}, {0}},
                   {-1, 2, {2, 1}, {0}},
                   {-1, 3, {3, 0}, {0}},
                   {-1, 3, {3, 1}, {0}}},
                  {// Tests for {0, 1, 0, 0, 0, 0, 0, 0} on the first rank,
                   // {0, 0, 2, 3, 0, 0, 0, 0} on the second, ...
                   {0, {1}},
                   {0, {0}},
                   {0, {0}},
                   {0, {4}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {0, {0}},
                   {1, {0}},
                   {1, {2}},
                   {1, {3}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {1, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {0}},
                   {2, {5}},
                   {2, {6}},
                   {2, {0}},
                   {2, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {0}},
                   {3, {7}},
                   {3, {8}}},
                  context.rank * 2);
}

/// Tests a non-contigous owner distributed at the outMesh
BOOST_AUTO_TEST_CASE(testDistributedConservative2DV5Vector)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  Gaussian                        fct(5.0);
  RadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  testDistributed(context, mapping,
                  {// Conservative mapping: The inMesh is local
                   {0, -1, {0, 0}, {1, 4}},
                   {0, -1, {0, 1}, {2, 5}},
                   {1, -1, {1, 0}, {3, 6}},
                   {1, -1, {1, 1}, {4, 7}},
                   {2, -1, {2, 0}, {5, 8}},
                   {2, -1, {2, 1}, {6, 9}},
                   {3, -1, {3, 0}, {7, 10}},
                   {3, -1, {3, 1}, {8, 11}}},
                  {// The outMesh is distributed and non-contigous
                   {-1, 0, {0, 0}, {0, 0}},
                   {-1, 1, {0, 1}, {0, 0}},
                   {-1, 1, {1, 0}, {0, 0}},
                   {-1, 0, {1, 1}, {0, 0}},
                   {-1, 2, {2, 0}, {0, 0}},
                   {-1, 2, {2, 1}, {0, 0}},
                   {-1, 3, {3, 0}, {0, 0}},
                   {-1, 3, {3, 1}, {0, 0}}},
                  {// Tests for {0, 1, 0, 0, 0, 0, 0, 0} on the first rank,
                   // {0, 0, 2, 3, 0, 0, 0, 0} on the second, ...
                   {0, {1, 4}},
                   {0, {0, 0}},
                   {0, {0, 0}},
                   {0, {4, 7}},
                   {0, {0, 0}},
                   {0, {0, 0}},
                   {0, {0, 0}},
                   {0, {0, 0}},
                   {1, {0, 0}},
                   {1, {2, 5}},
                   {1, {3, 6}},
                   {1, {0, 0}},
                   {1, {0, 0}},
                   {1, {0, 0}},
                   {1, {0, 0}},
                   {1, {0, 0}},
                   {2, {0, 0}},
                   {2, {0, 0}},
                   {2, {0, 0}},
                   {2, {0, 0}},
                   {2, {5, 8}},
                   {2, {6, 9}},
                   {2, {0, 0}},
                   {2, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {0, 0}},
                   {3, {7, 10}},
                   {3, {8, 11}}},
                  context.rank * 2);
}

void testTagging(const TestContext &context,
                 MeshSpecification  inMeshSpec,
                 MeshSpecification  outMeshSpec,
                 MeshSpecification  shouldTagFirstRound,
                 MeshSpecification  shouldTagSecondRound,
                 bool               consistent)
{
  int meshDimension  = inMeshSpec.at(0).position.size();
  int valueDimension = inMeshSpec.at(0).value.size();

  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", meshDimension, false, testing::nextMeshID()));
  mesh::PtrData inData = inMesh->createData("InData", valueDimension);
  getDistributedMesh(context, inMeshSpec, inMesh, inData);

  mesh::PtrMesh outMesh(new mesh::Mesh("outMesh", meshDimension, false, testing::nextMeshID()));
  mesh::PtrData outData = outMesh->createData("OutData", valueDimension);
  getDistributedMesh(context, outMeshSpec, outMesh, outData);

  Gaussian                        fct(4.5); //Support radius approx. 1
  Mapping::Constraint             constr = consistent ? Mapping::CONSISTENT : Mapping::CONSERVATIVE;
  RadialBasisFctMapping<Gaussian> mapping(constr, 2, fct, false, false, false);
  inMesh->computeBoundingBox();
  outMesh->computeBoundingBox();

  mapping.setMeshes(inMesh, outMesh);
  mapping.tagMeshFirstRound();

  for (const auto &v : inMesh->vertices()) {
    auto pos   = std::find_if(shouldTagFirstRound.begin(), shouldTagFirstRound.end(),
                            [meshDimension, &v](const VertexSpecification &spec) {
                              return std::equal(spec.position.data(), spec.position.data() + meshDimension, v.getCoords().data());
                            });
    bool found = pos != shouldTagFirstRound.end();
    BOOST_TEST(found >= v.isTagged(),
               "FirstRound: Vertex " << v << " is tagged, but should not be.");
    BOOST_TEST(found <= v.isTagged(),
               "FirstRound: Vertex " << v << " is not tagged, but should be.");
  }

  mapping.tagMeshSecondRound();

  for (const auto &v : inMesh->vertices()) {
    auto posFirst    = std::find_if(shouldTagFirstRound.begin(), shouldTagFirstRound.end(),
                                 [meshDimension, &v](const VertexSpecification &spec) {
                                   return std::equal(spec.position.data(), spec.position.data() + meshDimension, v.getCoords().data());
                                 });
    bool foundFirst  = posFirst != shouldTagFirstRound.end();
    auto posSecond   = std::find_if(shouldTagSecondRound.begin(), shouldTagSecondRound.end(),
                                  [meshDimension, &v](const VertexSpecification &spec) {
                                    return std::equal(spec.position.data(), spec.position.data() + meshDimension, v.getCoords().data());
                                  });
    bool foundSecond = posSecond != shouldTagSecondRound.end();
    BOOST_TEST(foundFirst <= v.isTagged(), "SecondRound: Vertex " << v
                                                                  << " is not tagged, but should be from the first round.");
    BOOST_TEST(foundSecond <= v.isTagged(), "SecondRound: Vertex " << v
                                                                   << " is not tagged, but should be.");
    BOOST_TEST((foundSecond or foundFirst) >= v.isTagged(), "SecondRound: Vertex " << v
                                                                                   << " is tagged, but should not be.");
  }
}

BOOST_AUTO_TEST_CASE(testTagFirstRound)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves())
  //    *
  //    + <-- owned
  //* * x * *
  //    *
  //    *
  MeshSpecification outMeshSpec = {
      {0, -1, {0, 0}, {0}}};
  MeshSpecification inMeshSpec = {
      {0, -1, {-1, 0}, {1}}, //inside
      {0, -1, {-2, 0}, {1}}, //outside
      {0, 0, {1, 0}, {1}},   //inside, owner
      {0, -1, {2, 0}, {1}},  //outside
      {0, -1, {0, -1}, {1}}, //inside
      {0, -1, {0, -2}, {1}}, //outside
      {0, -1, {0, 1}, {1}},  //inside
      {0, -1, {0, 2}, {1}}   //outside
  };
  MeshSpecification shouldTagFirstRound = {
      {0, -1, {-1, 0}, {1}},
      {0, -1, {1, 0}, {1}},
      {0, -1, {0, -1}, {1}},
      {0, -1, {0, 1}, {1}}};
  MeshSpecification shouldTagSecondRound = {
      {0, -1, {2, 0}, {1}}};
  testTagging(context, inMeshSpec, outMeshSpec, shouldTagFirstRound, shouldTagSecondRound, true);
  // For conservative just swap meshes.
  testTagging(context, outMeshSpec, inMeshSpec, shouldTagFirstRound, shouldTagSecondRound, false);
}

BOOST_AUTO_TEST_SUITE_END() // Parallel

BOOST_AUTO_TEST_SUITE(Serial)

void perform2DTestConsistentMapping(Mapping &mapping)
{
  int dimensions = 2;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector2d(0.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 1.0));
  inMesh->createVertex(Vector2d(0.0, 1.0));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  values << 1.0, 2.0, 2.0, 1.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector2d(0, 0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Vector2d(0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Vector2d(1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Vector2d(0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Vector2d(0.5, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Vector2d(0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);
}

void perform2DTestConsistentMappingVector(Mapping &mapping)
{
  int dimensions = 2;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 2);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector2d(0.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 0.0));
  inMesh->createVertex(Vector2d(1.0, 1.0));
  inMesh->createVertex(Vector2d(0.0, 1.0));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  values << 1.0, 4.0, 2.0, 5.0, 2.0, 5.0, 1.0, 4.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 2);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector2d(0, 0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Vector2d(0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value1 = outData->values()(0);
  double value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.0);
  BOOST_TEST(value2 == 4.0);

  vertex.setCoords(Vector2d(0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.0);
  BOOST_TEST(value2 == 4.0);

  vertex.setCoords(Vector2d(0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.0);
  BOOST_TEST(value2 == 4.0);

  vertex.setCoords(Vector2d(1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 2.0);
  BOOST_TEST(value2 == 5.0);

  vertex.setCoords(Vector2d(1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 2.0);
  BOOST_TEST(value2 == 5.0);

  vertex.setCoords(Vector2d(1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 2.0);
  BOOST_TEST(value2 == 5.0);

  vertex.setCoords(Vector2d(0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.5);
  BOOST_TEST(value2 == 4.5);

  vertex.setCoords(Vector2d(0.5, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.5);
  BOOST_TEST(value2 == 4.5);

  vertex.setCoords(Vector2d(0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value1 = outData->values()(0);
  value2 = outData->values()(1);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value1 == 1.5);
  BOOST_TEST(value2 == 4.5);
}

void perform3DTestConsistentMapping(Mapping &mapping)
{
  int dimensions = 3;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1);
  int           inDataID = inData->getID();
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  values << 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Eigen::Vector3d::Zero());
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(0.5, 0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(0.5, 0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value, 1.5);

  vertex.setCoords(Eigen::Vector3d(1.0, 0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(0.0, 1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(1.0, 1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(0.5, 0.5, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.5);
}

void perform2DTestConservativeMapping(Mapping &mapping)
{
  const int    dimensions = 2;
  const double tolerance  = 1e-6;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1);
  int           inDataID = inData->getID();
  mesh::Vertex &vertex0  = inMesh->createVertex(Vector2d(0, 0));
  mesh::Vertex &vertex1  = inMesh->createVertex(Vector2d(0, 0));
  inMesh->allocateDataValues();
  inData->values() << 1.0, 2.0;
  addGlobalIndex(inMesh);

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1);
  int           outDataID = outData->getID();
  outMesh->createVertex(Vector2d(0.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 1.0));
  outMesh->createVertex(Vector2d(0.0, 1.0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  auto &values = outData->values();

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex0.setCoords(Vector2d(0.5, 0.0));
  vertex1.setCoords(Vector2d(0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(testing::equals(values, Eigen::Vector4d(0.5, 0.5, 1.0, 1.0), tolerance));

  vertex0.setCoords(Vector2d(0.0, 0.5));
  vertex1.setCoords(Vector2d(1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(testing::equals(values, Eigen::Vector4d(0.5, 1.0, 1.0, 0.5), tolerance));

  vertex0.setCoords(Vector2d(0.0, 1.0));
  vertex1.setCoords(Vector2d(1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(testing::equals(values, Eigen::Vector4d(0.0, 2.0, 0.0, 1.0), tolerance));

  vertex0.setCoords(Vector2d(0.0, 0.0));
  vertex1.setCoords(Vector2d(1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(testing::equals(values, Eigen::Vector4d(1.0, 0.0, 2.0, 0.0), tolerance));

  vertex0.setCoords(Vector2d(0.4, 0.5));
  vertex1.setCoords(Vector2d(0.6, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(values.sum() == 3.0);
}

void perform2DTestConservativeMappingVector(Mapping &mapping)
{
  const int    dimensions = 2;
  const double tolerance  = 1e-6;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 2);
  int           inDataID = inData->getID();
  mesh::Vertex &vertex0  = inMesh->createVertex(Vector2d(0, 0));
  mesh::Vertex &vertex1  = inMesh->createVertex(Vector2d(0, 0));
  inMesh->allocateDataValues();
  inData->values() << 1.0, 4.0, 2.0, 5.0;
  addGlobalIndex(inMesh);

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 2);
  int           outDataID = outData->getID();
  outMesh->createVertex(Vector2d(0.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 0.0));
  outMesh->createVertex(Vector2d(1.0, 1.0));
  outMesh->createVertex(Vector2d(0.0, 1.0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  auto &values = outData->values();

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex0.setCoords(Vector2d(0.5, 0.0));
  vertex1.setCoords(Vector2d(0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  Eigen::VectorXd refValues(8);
  refValues << 0.5, 2, 0.5, 2, 1.0, 2.5, 1.0, 2.5;
  BOOST_TEST(testing::equals(values, refValues, tolerance));

  vertex0.setCoords(Vector2d(0.0, 0.5));
  vertex1.setCoords(Vector2d(1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  refValues << 0.5, 2, 1.0, 2.5, 1.0, 2.5, 0.5, 2;
  BOOST_TEST(testing::equals(values, refValues, tolerance));

  vertex0.setCoords(Vector2d(0.0, 1.0));
  vertex1.setCoords(Vector2d(1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  refValues << 0.0, 0.0, 2.0, 5.0, 0.0, 0.0, 1.0, 4.0;
  BOOST_TEST(testing::equals(values, refValues, tolerance));

  vertex0.setCoords(Vector2d(0.0, 0.0));
  vertex1.setCoords(Vector2d(1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  refValues << 1.0, 4.0, 0.0, 0.0, 2.0, 5.0, 0.0, 0.0;
  BOOST_TEST(testing::equals(values, refValues, tolerance));

  vertex0.setCoords(Vector2d(0.4, 0.5));
  vertex1.setCoords(Vector2d(0.6, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(values.sum() == 12.0);
}

void perform3DTestConservativeMapping(Mapping &mapping)
{
  using Eigen::Vector3d;
  int dimensions = 3;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1);
  int           inDataID = inData->getID();
  mesh::Vertex &vertex0  = inMesh->createVertex(Vector3d(0, 0, 0));
  mesh::Vertex &vertex1  = inMesh->createVertex(Vector3d(0, 0, 0));
  inMesh->allocateDataValues();
  inData->values() << 1.0, 2.0;
  addGlobalIndex(inMesh);

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1);
  int           outDataID = outData->getID();
  outMesh->createVertex(Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(0.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 1.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 1.0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  auto & values      = outData->values();
  double expectedSum = inData->values().sum();

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex0.setCoords(Vector3d(0.5, 0.0, 0.0));
  vertex1.setCoords(Vector3d(0.5, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping());
  BOOST_TEST(values.sum() == expectedSum);
}

BOOST_AUTO_TEST_CASE(MapThinPlateSplines)
{
  PRECICE_TEST(1_rank);
  bool                                    xDead = false;
  bool                                    yDead = false;
  bool                                    zDead = false;
  ThinPlateSplines                        fct;
  RadialBasisFctMapping<ThinPlateSplines> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<ThinPlateSplines> consistentMap2DVector(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMappingVector(consistentMap2DVector);
  RadialBasisFctMapping<ThinPlateSplines> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<ThinPlateSplines> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<ThinPlateSplines> conservativeMap2DVector(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMappingVector(conservativeMap2DVector);
  RadialBasisFctMapping<ThinPlateSplines> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapMultiquadrics)
{
  PRECICE_TEST(1_rank);
  bool                                 xDead = false;
  bool                                 yDead = false;
  bool                                 zDead = false;
  Multiquadrics                        fct(1e-3);
  RadialBasisFctMapping<Multiquadrics> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<Multiquadrics> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<Multiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<Multiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapInverseMultiquadrics)
{
  PRECICE_TEST(1_rank);
  bool                                        xDead = false;
  bool                                        yDead = false;
  bool                                        zDead = false;
  InverseMultiquadrics                        fct(1e-3);
  RadialBasisFctMapping<InverseMultiquadrics> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<InverseMultiquadrics> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<InverseMultiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<InverseMultiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapVolumeSplines)
{
  PRECICE_TEST(1_rank);
  bool                                 xDead = false;
  bool                                 yDead = false;
  bool                                 zDead = false;
  VolumeSplines                        fct;
  RadialBasisFctMapping<VolumeSplines> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<VolumeSplines> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<VolumeSplines> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<VolumeSplines> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapGaussian)
{
  PRECICE_TEST(1_rank);
  bool                            xDead = false;
  bool                            yDead = false;
  bool                            zDead = false;
  Gaussian                        fct(1.0);
  RadialBasisFctMapping<Gaussian> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<Gaussian> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<Gaussian> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<Gaussian> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapCompactThinPlateSplinesC2)
{
  PRECICE_TEST(1_rank);
  double                    supportRadius = 1.2;
  bool                      xDead         = false;
  bool                      yDead         = false;
  bool                      zDead         = false;
  CompactThinPlateSplinesC2 fct(supportRadius);
  using Mapping = RadialBasisFctMapping<CompactThinPlateSplinesC2>;
  Mapping consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapPetCompactPolynomialC0)
{
  PRECICE_TEST(1_rank);
  double              supportRadius = 1.2;
  bool                xDead         = false;
  bool                yDead         = false;
  bool                zDead         = false;
  CompactPolynomialC0 fct(supportRadius);
  using Mapping = RadialBasisFctMapping<CompactPolynomialC0>;
  Mapping consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapPetCompactPolynomialC6)
{
  PRECICE_TEST(1_rank);
  double              supportRadius = 1.2;
  bool                xDead         = false;
  bool                yDead         = false;
  bool                zDead         = false;
  CompactPolynomialC6 fct(supportRadius);
  using Mapping = RadialBasisFctMapping<CompactPolynomialC6>;
  Mapping consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(DeadAxis2)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector2d;
  int dimensions = 2;

  bool xDead = false;
  bool yDead = true;
  bool zDead = false;

  ThinPlateSplines                        fct;
  RadialBasisFctMapping<ThinPlateSplines> mapping(Mapping::CONSISTENT, dimensions, fct,
                                                  xDead, yDead, zDead);

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector2d(0.0, 1.0));
  inMesh->createVertex(Vector2d(1.0, 1.0));
  inMesh->createVertex(Vector2d(2.0, 1.0));
  inMesh->createVertex(Vector2d(3.0, 1.0));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  values << 1.0, 2.0, 2.0, 1.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1);
  int           outDataID = outData->getID();
  mesh::Vertex &vertex    = outMesh->createVertex(Vector2d(0, 0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex.setCoords(Vector2d(0.0, 3.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value = outData->values()(0);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(value == 1.0);
}

BOOST_AUTO_TEST_CASE(DeadAxis3D)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  int dimensions = 3;

  double              supportRadius = 1.2;
  CompactPolynomialC6 fct(supportRadius);
  bool                xDead = false;
  bool                yDead = true;
  bool                zDead = false;
  using Mapping             = RadialBasisFctMapping<CompactPolynomialC6>;
  Mapping mapping(Mapping::CONSISTENT, dimensions, fct, xDead, yDead, zDead);

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector3d(0.0, 3.0, 0.0));
  inMesh->createVertex(Vector3d(1.0, 3.0, 0.0));
  inMesh->createVertex(Vector3d(0.0, 3.0, 1.0));
  inMesh->createVertex(Vector3d(1.0, 3.0, 1.0));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);

  auto &values = inData->values();
  values << 1.0, 2.0, 3.0, 4.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1);
  int           outDataID = outData->getID();
  outMesh->createVertex(Vector3d(0.0, 2.9, 0.0));
  outMesh->createVertex(Vector3d(0.8, 2.9, 0.1));
  outMesh->createVertex(Vector3d(0.1, 2.9, 0.9));
  outMesh->createVertex(Vector3d(1.1, 2.9, 1.1));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outData->values()(0) == 1.0);
  BOOST_TEST(outData->values()(1) == 2.0);
  BOOST_TEST(outData->values()(2) == 2.9);
  BOOST_TEST(outData->values()(3) == 4.3);
}

BOOST_AUTO_TEST_SUITE_END() // Serial

BOOST_AUTO_TEST_SUITE_END() // RadialBasisFunctionMapping
BOOST_AUTO_TEST_SUITE_END()
