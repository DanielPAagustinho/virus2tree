import csv
import json
import argparse
from collections import defaultdict
from typing import Dict, Any
from typing import Tuple

def parse_csv_to_tree(input_file: str,tree_name: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    nodes = {}
    internal_node_count = 0
    leaf_node_count = 0
    children_map = defaultdict(list)
    # meta_keys = {}

    with open(input_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        types = next(reader)
        if not all(t in ["numeric", "integer", "Date","character"] for t in types.values()):
            raise ValueError("First row of CSV must contain type definitions like 'numeric', 'integer', or 'string'.")



        meta_keys = {k:{
            "type":v,
            "std_type":None,
            "name":k,
            "NA_count":0,
            "min":float("inf"),
            "max": float("-inf"),
            "categories":[]
        }  for k, v in types.items() if k not in ["parent", "node", "branch.length", "label","confidence"]}
        for row in reader:
            node_id = row["node"]
            parent_id = row["parent"]

            meta = {
                k: row[k]
                for k in row
                if k not in ["parent", "node", "branch.length", "label","confidence"]
            }

        

            for k, v in meta.items():
                if meta_keys[k]["type"] in ["numeric", "integer"]: # type: ignore
                    if v == "NA":
                        meta_keys[k]["NA_count"] += 1
                    else:
                        meta_keys[k]["min"] = min(meta_keys[k]["min"], float(v))
                        meta_keys[k]["max"] = max(meta_keys[k]["max"], float(v))
                    if meta_keys[k]["std_type"] is None:
                        meta_keys[k]["std_type"] = "cont"
                else:
                    if v == "NA":
                        meta_keys[k]["NA_count"] += 1
                    else:
                        if v not in meta_keys[k]["categories"]:
                            meta_keys[k]["categories"].append(v)
                    if meta_keys[k]["std_type"] is None:
                        meta_keys[k]["std_type"] = "cate"

            brlen = 0.0 if row["branch.length"] == "NA" else float(row["branch.length"])

            nodes[node_id] = {
                "id": node_id,
                "name": row["label"],
                "branch_length": brlen,
                "meta": meta,
                "confidence": row["confidence"],
                "children": []
            }

            if parent_id != node_id:
                children_map[parent_id].append(node_id)
            else:
                root_id = node_id
        
            # meta_keys = list(meta.keys())
        # print("Meta keys:", json.dumps(meta_keys, indent=2))
        


    def build_tree(node_id: str) -> Dict[str, Any]:
        node = nodes[node_id]
        node["children"] = [build_tree(child_id) for child_id in children_map.get(node_id, [])]
        return node
    
    if 'root_id' not in locals():
        raise ValueError("Root node not found. Ensure there is a row where 'parent' == 'node' to define the root.")

    tree = build_tree(root_id)


    for node_id, node in nodes.items():
        if len(node['children']) == 0:
            leaf_node_count += 1
        else:
            internal_node_count += 1

    tree_meta = {
        "name": tree_name,
        "node_count": len(nodes),
        "leaf_count": leaf_node_count,
        "internal_count": internal_node_count,
        "meta_keys": meta_keys
    }

    return  (tree, tree_meta)

def estimate_tree_dimensions(tree: Dict[str, Any]) -> Dict[str, int]:
    def dfs(node, depth):
        nonlocal max_depth, leaf_count
        max_depth = max(max_depth, depth)
        if not node.get("children"):  # 叶子节点
            leaf_count += 1
        for child in node.get("children", []):
            dfs(child, depth + 1)

    max_depth = 0
    leaf_count = 0
    dfs(tree, 1)

    return {
        "estimated_width": max_depth,
        "estimated_height": leaf_count
    }

def main(input_file: str, out_prefix: str, html_template: str, tree_name: str):
    tree,tree_meta = parse_csv_to_tree(input_file, tree_name)


    # clean the inf and -inf values

    for k, v in tree_meta["meta_keys"].items():
        if v["min"] == float("inf"):
            v["min"] = None
        if v["max"] == float("-inf"):
            v["max"] = None

    # estimate the tree dimensions
    dimensions = estimate_tree_dimensions(tree)
    tree_meta["estimated_width"] = dimensions["estimated_width"]
    tree_meta["estimated_height"] = dimensions["estimated_height"]


    out_html = out_prefix + ".html"
    tree_json_output = out_prefix + ".tree.json"
    tree_meta_json_output = out_prefix + ".tree_meta.json"
    with open(tree_json_output, "w") as f:
        json.dump(tree, f, indent=2)
    with open(tree_meta_json_output, "w") as f2:
        json.dump(tree_meta, f2, indent=2)
    tree_str = json.dumps(tree)
    with open(html_template, "r") as f:
        template = f.read()
    with open(out_html, "w") as f2:

        template = template.replace("[/*inject_meta*/]", json.dumps(tree_meta))

        f2.write(template.replace("[/*inject_data*/]", tree_str))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert phylogenetic CSV to JSON tree and inject into HTML.")
    parser.add_argument("-i","--input_file",required=True, help="Path to input CSV file")
    parser.add_argument("-o","--out",required=True, help="Path of output prefix")
    parser.add_argument("-t","--template",required=True, help="Path to input HTML template")
    parser.add_argument("-n","--name",required=True, help="Name of the tree",default="tree")
    
    args = parser.parse_args()
    if not args.input_file or not args.out or not args.template:
        parser.print_help()
        exit(1)
    main(args.input_file, args.out, args.template, args.name)
