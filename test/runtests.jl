using Test
using Random

using EpiSim
using Outbreak: Outbreak,
                has_log,
                has_tree,
                has_alignment,
                iscomplete,
                with_tree,
                with_alignment,
                outbreak_tree,
                outbreak_alignment,
                simulate_outbreak_tree,
                simulate_outbreak_alignment
using SeqSim
using TreeSim

function fixture_event_log()
    return EventLog(
        [0.0, 1.0, 2.0, 3.0, 4.0],
        [1, 2, 3, 2, 3],
        [0, 1, 1, 0, 0],
        [EK_Seeding, EK_Transmission, EK_Transmission, EK_SerialSampling, EK_SerialSampling],
    )
end

function fixture_tree()
    return Tree(
        [0.0, 0.4, 0.8, 1.1, 1.3],
        [2, 4, 0, 0, 0],
        [3, 5, 0, 0, 0],
        [0, 1, 1, 2, 2],
        [Root, Binary, SampledLeaf, SampledLeaf, SampledLeaf],
        [0, 0, 0, 0, 0],
        [900, 0, 103, 104, 105],
    )
end

site_model() = SiteModel(MersenneTwister(11), 12, 0.2, 0, 0.0, 0.0, JC())

function same_tree(a::Tree, b::Tree)
    return a.time == b.time &&
           a.left == b.left &&
           a.right == b.right &&
           a.parent == b.parent &&
           a.kind == b.kind &&
           a.host == b.host &&
           a.label == b.label
end

@testset "Outbreak bundle construction" begin
    log = fixture_event_log()
    tree = fixture_tree()
    aln = [Sequence("ACGT"; taxon="a"), Sequence("ACGA"; taxon="b")]

    from_log = Outbreak(log)
    from_tree = Outbreak(log, tree)
    full = Outbreak(log, tree, aln)

    @test from_log.log === log
    @test from_log.tree === nothing
    @test from_log.aln === nothing
    @test has_log(from_log)
    @test !has_tree(from_log)
    @test !has_alignment(from_log)
    @test !iscomplete(from_log)

    @test from_tree.tree === tree
    @test from_tree.aln === nothing
    @test has_tree(from_tree)
    @test !has_alignment(from_tree)
    @test !iscomplete(from_tree)

    @test full.aln === aln
    @test has_alignment(full)
    @test iscomplete(full)
    @test sprint(show, full) == "Outbreak(log=EventLog, tree=Tree, aln=Vector)"
    @test occursin("log  : populated", sprint(show, MIME"text/plain"(), full))
end

@testset "Event log to tree orchestration" begin
    log = fixture_event_log()

    out = with_tree(log)
    simulated = simulate_outbreak_tree(log)

    @test out.log === log
    @test has_tree(out)
    @test same_tree(out.tree, tree_from_eventlog(log))
    @test validate_tree(out.tree)
    @test validate_tree_against_eventlog(log, out.tree)
    @test same_tree(outbreak_tree(log), out.tree)
    @test simulated.log === log
    @test same_tree(simulated.tree, out.tree)
    @test simulated.aln === nothing
end

@testset "Invalid stage transitions fail explicitly" begin
    log = fixture_event_log()
    tree = fixture_tree()
    model = site_model()

    tree_only = Outbreak(nothing, tree)
    log_only = Outbreak(log)
    empty_tree = Outbreak(nothing, nothing)

    err = @test_throws ArgumentError with_tree(tree_only)
    @test occursin("log stage is empty", sprint(showerror, err.value))

    err = @test_throws ArgumentError with_tree(empty_tree)
    @test occursin("populate the tree stage", sprint(showerror, err.value))

    err = @test_throws ArgumentError outbreak_tree(nothing)
    @test occursin("log stage is empty", sprint(showerror, err.value))

    err = @test_throws ArgumentError with_alignment(log_only, model)
    @test occursin("tree stage is empty", sprint(showerror, err.value))

    err = @test_throws ArgumentError with_alignment(MersenneTwister(91), empty_tree, model)
    @test occursin("populate the alignment stage", sprint(showerror, err.value))
end

@testset "Tree to alignment orchestration" begin
    tree = fixture_tree()
    model = site_model()

    out = with_alignment(MersenneTwister(21), tree, model)
    direct = simulate_alignment(MersenneTwister(21), tree, model)

    @test out.log === nothing
    @test out.tree === tree
    @test out.aln == direct
    @test length(out.aln) == nleaves(tree)
    @test outbreak_alignment(MersenneTwister(21), tree, model) == direct
    @test_throws ArgumentError with_alignment(Outbreak(fixture_event_log()), model)

    @test sprint(show, out) == "Outbreak(log=empty, tree=Tree, aln=Vector)"
    @test occursin("log  : empty", sprint(show, MIME"text/plain"(), out))
end

@testset "End-to-end smoke path" begin
    log = fixture_event_log()
    model = site_model()

    out = simulate_outbreak_alignment(MersenneTwister(31), log, model)
    expected_tree = tree_from_eventlog(log)
    expected_aln = simulate_alignment(MersenneTwister(31), expected_tree, model)

    @test out.log === log
    @test same_tree(out.tree, expected_tree)
    @test out.aln == expected_aln
    @test length(out.aln) == nleaves(out.tree)
    @test [seq.taxon for seq in out.aln] == out.tree.label[tips(out.tree)]
end

@testset "Default RNG overloads use task-local default RNG" begin
    log = fixture_event_log()
    tree = fixture_tree()
    model = site_model()
    rng = Random.default_rng()

    Random.seed!(rng, 71)
    default_tree_out = with_alignment(tree, model)
    Random.seed!(rng, 71)
    explicit_tree_out = with_alignment(rng, tree, model)
    @test default_tree_out.aln == explicit_tree_out.aln

    Random.seed!(rng, 72)
    default_aln = outbreak_alignment(tree, model)
    Random.seed!(rng, 72)
    explicit_aln = outbreak_alignment(rng, tree, model)
    @test default_aln == explicit_aln

    Random.seed!(rng, 73)
    default_full = simulate_outbreak_alignment(log, model)
    Random.seed!(rng, 73)
    explicit_full = simulate_outbreak_alignment(rng, log, model)
    @test default_full.aln == explicit_full.aln
    @test same_tree(default_full.tree, explicit_full.tree)

    @test !occursin("GLOBAL_RNG", read(joinpath(dirname(@__DIR__), "src", "Outbreak.jl"), String))
end
