#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use threads;

# Copyright (c) 2021 Hikoyu Suzuki
# This software is released under the MIT License.

# ソフトウェアを定義
### 編集範囲 開始 ###
my $software = "divgraph.pl";	# ソフトウェアの名前
my $version = "ver.1.0.0";	# ソフトウェアのバージョン
my $note = "Divide a weighted graph into communities by using greedy method under specified conditions.";	# ソフトウェアの説明
my $usage = "<STDIN|edge1.tsv> [edge2.tsv ...] [>cluster.tsv]";	# ソフトウェアの使用法 (コマンド非使用ソフトウェアの時に有効)
### 編集範囲 終了 ###

# コマンドを定義
my %command;
### 編集範囲 開始 ###
# コマンドを追加
### 編集範囲 終了 ###
my @command_list = sort(keys(%command));

# 指定されたコマンドを確認
my $specified_command = shift(@ARGV) if @command_list and @ARGV;
&exception::error("unknown command: $specified_command") if $specified_command and !grep {$_ eq $specified_command} @command_list;

# 共通オプションを定義
my %option;
### 編集範囲 開始 ###
$option{"p INT "} = "Number of parallel worker threads <1-> [1]";
$option{"t INT "} = "Number of trials per worker thread for clustering sequence reads <1-> [1]";
$option{"w"} = "\tUse 2-byte line feed code (CR+LF) for input files";
# オプションを追加
### 編集範囲 終了 ###

# コマンドごとのオプション定義を取得
&{\&{"${specified_command}::define"}} if $specified_command;
my @option_list = sort(keys(%option));

# ヘルプを表示 (引数未指定時)
&exception::help if !@ARGV and !-p STDIN;

# オプションの入力処理
my %opt;
$_ = join("", @option_list);
$_ =~ s/\s+\S+\s+/:/g;
getopts($_, \%opt);

# 未指定オプションのデフォルト値を入力
foreach (@option_list) {
	$opt{substr($_, 0, 1)} = substr($option{$_}, index($option{$_}, "[") + 1, index($option{$_}, "]") - index($option{$_}, "[") - 1) if $option{$_} =~ /\[.+\]$/ and !defined($opt{substr($_, 0, 1)});
}

### 編集範囲 開始 ###
# 追加のモジュールを宣言
use threads::shared;
use Thread::Queue 3.07;
use Inline (CPP => Config => CC => exists($ENV{"CXX"}) ? $ENV{"CXX"} : 'c++', CCFLAGS => '-std=c++14 -march=native', DIRECTORY => $ENV{"DIVGRAPH_INLINE_DIR"});
use Inline (CPP => 'DATA', NAME => 'DIVGRAPH', AUTO_INCLUDE => ['#undef seed', '#include <vector>', '#include <unordered_map>', '#include <random>', '#include <numeric>']);
no warnings 'portable';

# 定数を定義
use constant max_nodes => 4294967295;	# max_nodes => ノード数上限値

# 処理を追加
### 編集範囲 終了 ###

# メインルーチンを実行
&main;
exit(0);

# メインルーチン
sub main {
	### 編集範囲 開始 ###
	# 指定されたオプションを確認
	&exception::error("specify INT >= 1: -p $opt{p}") if $opt{"p"} !~ /^\d+$/ or $opt{"p"} < 1;
	&exception::error("specify INT >= 1: -t $opt{t}") if $opt{"t"} !~ /^\d+$/ or $opt{"t"} < 1;
	
	# 変数を宣言
	my %node_list = ();
	my @node_name = ();
	my @edge_list = ();
	
	# 入力の改行コードを一時的に変更 (-w指定時)
	local $/ = "\r\n" if $opt{"w"};
	
	# エッジリストデータを読み込みながら処理
	print STDERR "Loading edge lists...";
	while (<>) {
		# 改行コードを除去
		chomp;
		
		# 空白文字で分割
		my ($node1, $node2, $weight) = split(/\s+/);
		
		# ノードが未定義の場合は定義する
		push(@node_name, $node1) and $node_list{$node1} = @node_name unless $node_list{$node1};
		push(@node_name, $node2) and $node_list{$node2} = @node_name unless $node_list{$node2};
		
		# 重みが未定義の場合は1とする
		$weight = 1 unless defined($weight);
		
		# 重みが0でない場合はエッジリストに登録
		$edge_list[$node_list{$node1}] .= pack("Lf", $node_list{$node2}, $weight) and $edge_list[$node_list{$node2}] .= pack("Lf", $node_list{$node1}, $weight) if $weight;
	}
	print STDERR "completed\n";
	
	# ノードが存在しない場合は除外
	&exception::error("No nodes in the input graph data.") unless %node_list;
	
	# ノード数が上限値を超える場合は除外
	&exception::error("Too many nodes in the input graph data.") if @node_name > max_nodes;
	
	# 変数を宣言
	my @max_clusters = ();
	my $already_cluster_assigned = "";
	
	### 互いに到達可能なノードのみからなる極大クラスターを列挙 ###
	print STDERR "Creating maximum clusters...";
	# 各ノードについて到達可能なノードを列挙し極大クラスターを決定
	for (my $node_id = 1;$node_id <= @node_name;$node_id++) {
		# クラスターに割り当て済みの場合は以下の処理をスキップ
		next if vec($already_cluster_assigned, $node_id, 1);
		
		# 変数を宣言
		my @node_queue = ($node_id);
		
		# クラスターを作成
		my $max_cluster = pack("L", $node_id);
		
		# クラスター割り当てフラグをチェック
		vec($already_cluster_assigned, $node_id, 1) = 1;
		
		# ノードキューからノードIDを取得して処理
		while (my $node_id = shift(@node_queue)) {
			# ノードIDを取得
			my @node_ids = unpack("(Lx4)*", $edge_list[$node_id]);
			
			# クラスターに未割り当てのノードについてノードキューとクラスターに追加してクラスター割り当てフラグをチェック
			push(@node_queue, $_) and $max_cluster .= pack("L", $_) and vec($already_cluster_assigned, $_, 1) = 1 foreach grep {!vec($already_cluster_assigned, $_, 1)} @node_ids;
		}
		
		# 作成したクラスターをバイナリ形式で極大クラスターリストに追加
		push(@max_clusters, $max_cluster);
	}
	print STDERR "completed\n";
	
	### excess modularity density (Qx) に基づくクラスタリング ###
	# 共有変数を宣言
	my $max_Qx : shared;
	my $best_cluster_assignment : shared;
	
	# Qx最大値の初期値を設定
	$max_Qx = -Inf;
	
	# 変数を宣言
	my @worker_threads = ();
	my $num_error_threads = 0;
	
	# 指定されたワーカースレッド数で並列処理
	print STDERR "Creating optimized clusters...";
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		$worker_threads[$thread_id] = threads::async {
			# このスレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 試行回数だけ処理
			for (my $i = 0;$i < $opt{"t"};$i++) {
				# クラスター割り当てリストを初期化
				my $cluster_assignment = pack("L*", 0..@node_name);
				
				# excess modularity density (Qx) に基づいてノードをクラスタリング
				my $Qx = &create_clusters($cluster_assignment, \@max_clusters, \@edge_list);
				
				# 共有変数をロック
				lock($max_Qx);
				
				# Qxが最大値を超えた場合はQx最大値及び最善のクラスター割り当てリストを更新
				($max_Qx, $best_cluster_assignment) = ($Qx, $cluster_assignment) if $Qx > $max_Qx;
			}
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# 各ワーカースレッドが終了するまで待機
	$worker_threads[$_]->join or $num_error_threads++ foreach 0..$opt{"p"} - 1;
	
	# 異常終了したワーカースレッド数を確認
	print STDERR "aborted\n" and &exception::error("$num_error_threads worker thread" . ($num_error_threads > 1 ? "s" : "") . " abnormally exited") if $num_error_threads;
	print STDERR "completed\n";
	
	# モジュラリティを出力
	print STDERR "Best excess modularity density: $max_Qx\n";
	
	# 変数を宣言
	my %members = ();
	my $num_clusters = 0;
	
	# 各ノードについて所属クラスターを登録
	for (my $i = 1;$i <= @node_name;$i++) {push(@{$members{substr($best_cluster_assignment, $i * 4, 4)}}, $i - 1);}
	
	# 各クラスターについて所属ノードを出力
	print join("\t", "Node ID", "Cluster ID"), "\n";
	++$num_clusters and print join("\n", map {join("\t", $node_name[$_], $num_clusters)} @{$members{$_}}), "\n" foreach sort {@{$members{$b}} <=> @{$members{$a}} || $a cmp $b} keys(%members);
	### 編集範囲 終了 ###
	
	return(1);
}

## ここから例外処理のパッケージ ##
package exception;

# ヘルプ表示
sub help {
	print STDERR "$software ";
	print STDERR $specified_command ? $specified_command : $version;
	print STDERR "\n\nFunctions:\n  $note\n\nUsage:\n  $software ";
	if (!$specified_command and @command_list) {
		print STDERR "<command>\n";
		print STDERR "\nCommand:\n";
		foreach (@command_list) {print STDERR "  $_\t$command{$_}\n";}
	}
	else {
		print STDERR "$specified_command " if $specified_command;
		print STDERR "[options] " if @option_list;
		print STDERR "$usage\n";
		print STDERR "\nOptions:\n" if @option_list;
		foreach (@option_list) {print STDERR "  -$_\t$option{$_}\n";}
	}
	exit(0);
}

# エラー表示
sub error {
	print STDERR $software;
	print STDERR " $specified_command" if $specified_command;
	print STDERR ": Error: $_[0]";
	print STDERR ": $_[1] line $." if $_[1];
	print STDERR "\n";
	#threads->tid or map {$_->detach} threads->list;	# threadsモジュールを使用する場合はアンコメント
	exit(1);
}

# 注意表示
sub caution {
	print STDERR $software;
	print STDERR " $specified_command" if $specified_command;
	print STDERR ": Caution: $_[0]";
	print STDERR ": $_[1] line $." if $_[1];
	print STDERR "\n";
	return(1);
}

### 編集範囲 開始 ###
# サブルーチンを追加

# パッケージを追加
### 編集範囲 終了 ###
__END__
__CPP__
// 構造体を宣言
typedef struct _edge {
	uint32_t destination;
	float weight;
} edge;

// 密度を算出 main::calc_density(内部エッジの重みの総和, ノード数)
float calc_density(float sum_internal_edges, unsigned int num_nodes) {
	// 密度を算出
	float density = num_nodes > 1 ? sum_internal_edges / num_nodes / (num_nodes - 1) : 0.0f;
	
	// 密度を返す
	return density;
}

// 部分的なQxを算出 main::calc_partial_Qx(内部エッジの重みの総和, ノード数, クラスター次数, グラフ次数, グラフ密度)
float calc_partial_Qx(float sum_internal_edges, unsigned int num_nodes, float cluster_degree, float graph_degree, float graph_density) {
	// 密度指数を算出
	float density_index = calc_density(sum_internal_edges, num_nodes) - graph_density;
	
	// 部分的なQxを算出
	float partial_Qx = sum_internal_edges * density_index - pow(cluster_degree * density_index, 2) / graph_degree;
	
	// 部分的なQxを返す
	return partial_Qx;
}

// excess modularity density (Qx) に基づくノードのクラスタリング main::create_clusters(クラスター割り当てリストリファレンス, 極大クラスターリストリファレンス, エッジリストリファレンス)
float create_clusters(SV* cluster_assignment_sv, AV* max_clusters_av, AV* edge_list_av) {
	// クラスター割り当てリストを取得
	uint32_t* cluster_assignment = (uint32_t*)SvPV_nolen(cluster_assignment_sv);
	
	// ノード数を取得
	uint32_t num_nodes = SvCUR(cluster_assignment_sv) / sizeof(uint32_t) - 1;
	
	// 変数を宣言
	std::vector<std::vector<uint32_t>> max_clusters(av_len(max_clusters_av) + 1);
	std::vector<edge*> edge_list(num_nodes + 1);
	std::vector<uint32_t> num_edges(num_nodes + 1);
	std::vector<float> node_degree(num_nodes + 1, 0.0f);
	std::vector<float> sum_internal_edges(num_nodes + 1, 0.0f);
	std::vector<uint32_t> cluster_size(num_nodes + 1, 1);
	std::vector<float> cluster_degree(num_nodes + 1, 0.0f);
	std::vector<uint32_t> node_queue;
	std::vector<bool> cluster_updated(num_nodes + 1, 0);
	float last_Qx = log(0.0f);
	
	// 各極大クラスターについて処理
	for (uint32_t i = 0;i < max_clusters.size();i++) {
		SV* max_cluster_sv = *av_fetch(max_clusters_av, i, 0);
		uint32_t* max_cluster = (uint32_t*)SvPV_nolen(max_cluster_sv);
		max_clusters[i].assign(max_cluster, max_cluster + SvCUR(max_cluster_sv) / sizeof(uint32_t));
	}
	
	// 変数を初期化
	for (uint32_t node_id = 1;node_id <= num_nodes;node_id++) {
		// エッジリストを取得
		SV* edge_list_sv = *av_fetch(edge_list_av, node_id, 0);
		edge_list[node_id] = (edge*)SvPV_nolen(edge_list_sv);
		
		// エッジ数を取得
		num_edges[node_id] = SvCUR(edge_list_sv) / sizeof(edge);
		
		// 重みの和を取得
		node_degree[node_id] = std::accumulate(edge_list[node_id], edge_list[node_id] + num_edges[node_id], 0.0f, [](float x, edge y) {return x + y.weight;});
		cluster_degree[node_id] = node_degree[node_id];
		
		// ノードキューを初期化
		node_queue.push_back(node_id);
	}
	
	// グラフ次数を算出
	float graph_degree = std::accumulate(cluster_degree.begin(), cluster_degree.end(), 0.0f);
	
	// グラフ密度を算出
	float graph_density = calc_density(graph_degree, num_nodes);
	
	// 初期Qxを算出
	float Qx = std::accumulate(cluster_degree.begin(), cluster_degree.end(), 0.0f, [graph_degree, graph_density](float x, float y) {return x + calc_partial_Qx(0.0f, 1, y, graph_degree, graph_density);});
	
	// メルセンヌ・ツイスターによる擬似乱数生成器を作成
	std::random_device seed;
	std::mt19937 mt(seed());
	
	// Qxが向上した場合は処理を継続
	while (Qx > last_Qx) {
		// 現在のQxを保存
		last_Qx = Qx;
		
		// クラスター更新フラグをリセット
		cluster_updated.assign(cluster_updated.size(), 0);
		
		// 各ノードについて処理
		for (uint32_t i = node_queue.size() - 1;i < node_queue.size();i--) {
			// 0以上i以下の整数値の一様分布器を作成
			std::uniform_int_distribution<> rnd(0, i);
			
			// ランダムにノードを選択
			uint32_t j = rnd(mt);
			uint32_t node = node_queue[j];
			node_queue[j] = node_queue[i];
			node_queue[i] = node;
			
			// 所属クラスターを取得
			uint32_t cluster = cluster_assignment[node];
			
			// 変数を宣言
			std::unordered_map<uint32_t, float> sum_linked_edges{{cluster, 0.0f}, {0, 0.0f}};
			
			// 隣接クラスターに対するエッジの重みの総和を算出
			for (uint32_t j = 0;j < num_edges[node];j++) {sum_linked_edges[cluster_assignment[edge_list[node][j].destination]] += 2.0f * edge_list[node][j].weight;}
			
			// 所属クラスターから転出する場合のQx変化量を算出
			float basal_Qx_change = calc_partial_Qx(sum_internal_edges[cluster] - sum_linked_edges[cluster], cluster_size[cluster] - 1, cluster_degree[cluster] - node_degree[node], graph_degree, graph_density) - calc_partial_Qx(sum_internal_edges[cluster], cluster_size[cluster], cluster_degree[cluster], graph_degree, graph_density);
			
			// 隣接クラスターに転入 (あるいは所属クラスターから独立) する場合のQx変化量を算出
			std::unordered_map<uint32_t, float> Qx_change;
			for (const auto& sum_linked_edge : sum_linked_edges) {
				Qx_change[sum_linked_edge.first] = calc_partial_Qx(sum_internal_edges[sum_linked_edge.first] + sum_linked_edge.second, cluster_size[sum_linked_edge.first] + 1, cluster_degree[sum_linked_edge.first] + node_degree[node], graph_degree, graph_density) - calc_partial_Qx(sum_internal_edges[sum_linked_edge.first], cluster_size[sum_linked_edge.first], cluster_degree[sum_linked_edge.first], graph_degree, graph_density);
			}
			
			// 所属クラスターを変えない場合のQx変化量を定義
			Qx_change[cluster] = -basal_Qx_change;
			
			// Qx変化量が最大となるクラスターを選出
			uint32_t target_cluster = std::max_element(Qx_change.begin(), Qx_change.end(), [](const auto& x, const auto& y) {return x.second < y.second | x.second == y.second & x.first < y.first;})->first;
			
			// 所属クラスターから独立する場合は新規クラスターを登録
			if (!target_cluster) {
				while (cluster_size[target_cluster]) {target_cluster++;}
				sum_linked_edges[target_cluster] = 0.0f;
				Qx_change[target_cluster] = Qx_change[0];
			}
			
			// クラスター内エッジの重みの総和を更新
			sum_internal_edges[target_cluster] += sum_linked_edges[target_cluster];
			sum_internal_edges[cluster] -= sum_linked_edges[cluster];
			
			// クラスターサイズを更新
			cluster_size[target_cluster]++;
			cluster_size[cluster]--;
			
			// クラスター次数を更新
			cluster_degree[target_cluster] += node_degree[node];
			cluster_degree[cluster] -= node_degree[node];
			
			// 所属クラスターを更新
			cluster_assignment[node] = target_cluster;
			
			// Qxを更新
			Qx += Qx_change[target_cluster] + basal_Qx_change;
			
			// クラスター更新フラグをチェック
			cluster_updated[target_cluster] = target_cluster != cluster;
			cluster_updated[cluster] = target_cluster != cluster;
		}
		
		// 各極大クラスターについて属するノードの所属クラスターのいずれかが更新されたか否かで分類
		auto part = std::partition(max_clusters.begin(), max_clusters.end(), [cluster_updated, cluster_assignment](const auto& x) {return std::any_of(x.begin(), x.end(), [=](uint32_t y) {return cluster_updated[cluster_assignment[y]];});});
		
		// 所属クラスターが固定されたノードを取得
		std::unordered_map<uint32_t, uint32_t> fixed_nodes;
		for (auto max_cluster = part;max_cluster != max_clusters.end();max_cluster++) {
			for (uint32_t i = 0;i < (*max_cluster).size();i++) {
				fixed_nodes[(*max_cluster)[i]] = 0;
			}
		}
		
		// 属するノードの所属クラスターがいずれも固定された極大クラスターを削除
		max_clusters.erase(part, max_clusters.end());
		
		// 所属クラスターが固定された各ノードについてノードキューの位置を取得
		for (uint32_t i = 0;i < node_queue.size();i++) {
			if (fixed_nodes.find(node_queue[i]) != fixed_nodes.end()) {fixed_nodes[node_queue[i]] = i;}
		}
		
		// 固定されたノード数を取得
		uint32_t num_fixed_nodes = fixed_nodes.size();
		
		// 変数を宣言
		auto fixed_node = fixed_nodes.begin();
		
		// 所属クラスターが固定された各ノードについてノードキューの末尾に配置
		for (uint32_t i = node_queue.size() - num_fixed_nodes;i < node_queue.size();i++, fixed_node++) {
			std::swap(node_queue[fixed_node->second], node_queue[i]);
			if (fixed_nodes.find(node_queue[fixed_node->second]) != fixed_nodes.end()) {fixed_nodes[node_queue[fixed_node->second]] = fixed_node->second;}
		}
		
		// 所属クラスターが固定された各ノードをノードキューから削除
		node_queue.erase(node_queue.end() - num_fixed_nodes, node_queue.end());
	}
	
	// Qxを返す
	return Qx / graph_degree;
}
