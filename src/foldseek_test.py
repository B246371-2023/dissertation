# import os
# import sys
# import time
# import pandas as pd
# import json
# from requests import get, post
import os
import sys
import time
import pandas as pd
import json
import asyncio
import aiohttp
from requests import get, post,RequestException
import subprocess

# 示例PDB文件路径
pdb_path = "/home/s2530615/fold_tree/test/structs/A0A2J8ENW0.pdb"
output_file = "output.json"

def foldseek_search(pdb_path, output_file):
    try:
        # Load PDB
        time.sleep(5)
        if not os.path.exists(output_file):
            print("Reading PDB file...")
            with open(pdb_path, 'r') as file:
                data = file.read()
            print("PDB file read successfully. Starting search...")

            # Starting search
            ticket = post('https://search.foldseek.com/api/ticket', {
                            'q': data,
                            'database[]': ['afdb50'],
                            'mode': '3diaa',
                        }).json()
            print("Search started. Ticket ID:", ticket['id'])

            # Collecting results
            repeat = True
            while repeat:
                status = get('https://search.foldseek.com/api/ticket/' + ticket['id']).json()
                if status['status'] == "ERROR":
                    # Handle error
                    print("Error in status:", status)
                    sys.exit(0)
                print("Status:", status)

                # Wait a short time between poll requests
                time.sleep(1)
                repeat = status['status'] != "COMPLETE"
            print("Search complete. Retrieving results...")

            # Get all hits for the first query (0)
            #result = get('https://search.foldseek.com/api/result/' + ticket['id'] + '/0',timeout=60)
            subprocess.run(['wget', '-O', 'data.json1', f'https://search.foldseek.com/api/result/{ticket["id"]}/0'])
            print("Results retrieved successfully.")
            
            
            result.raise_for_status()  # 检查 HTTP 请求是否成功
            print("Results retrieved successfully.")
            data = result.json()
                
            with open('data.json', 'w') as file:
                json.dump(data, file, indent=4)
            
            print("JSON 数据已成功下载并保存到 data.json 文件。")

            # Iterating over results
            foldseek_results_rows = []
            for results_dbs in result['results']:
                for alignment in results_dbs['alignments']:
                    foldseek_results_rows.append(pd.DataFrame.from_dict({'pdb_file': [pdb_path],
                                                                         'db': [results_dbs['db']], 
                                                                         'query': [alignment.get('query', '-')],
                                                                         'target': [alignment.get('target', '-')],
                                                                         'seqId': [alignment.get('seqId', '-')],
                                                                         'alnLength': [alignment.get('alnLength', '-')],
                                                                         'missmatches': [alignment.get('missmatches', '-')],
                                                                         'gapsopened': [alignment.get('gapsopened', '-')],
                                                                         'qStartPos': [alignment.get('qStartPos', '-')],
                                                                         'qEndPos': [alignment.get('qEndPos', '-')],
                                                                         'dbStartPos': [alignment.get('dbStartPos', '-')],
                                                                         'dbEndPos': [alignment.get('dbEndPos', '-')],
                                                                         'prob': [alignment.get('prob', '-')],
                                                                         'eval': [alignment.get('eval', '-')],
                                                                         'score': [alignment.get('score', '-')],
                                                                         'qLen': [alignment.get('qLen', '-')],
                                                                         'dbLen': [alignment.get('dbLen', '-')],
                                                                         'qAln': [alignment.get('qAln', '-')],
                                                                         'dbAln': [alignment.get('dbAln', '-')],
                                                                         'tCa': [alignment.get('tCa', '-')],
                                                                         'tSeq': [alignment.get('tSeq', '-')],
                                                                         'taxId': [alignment.get('taxId', '-')],
                                                                         'taxName': [alignment.get('taxName', '-')]}))
            return pd.concat(foldseek_results_rows)
        else:
            print("Output file already exists.")
    except Exception as e:
        print(e)

# 测试函数
result_df = foldseek_search(pdb_path, output_file)
if result_df is not None:
    print(result_df.head())
else:
    print("No results returned.")

# async def fetch(session, url):
#     async with session.get(url) as response:
#         response_text = await response.text()
#         print(f"Response from {url}: {response_text[:500]}")  # 打印前500个字符进行调试
#         if response.status != 200:
#             raise ClientResponseError(response.request_info, response.history, status=response.status)
#         try:
#             return await response.json()
#         except aiohttp.ContentTypeError:
#             print("Response is not in JSON format")
#             return await response_text

# async def foldseek_search(pdb_path, output_file, retries=5, timeout=60):
#     try:
#         # Load PDB
#         await asyncio.sleep(5)
#         if not os.path.exists(output_file):
#             print("Reading PDB file...")
#             with open(pdb_path, 'r') as file:
#                 data = file.read()
#             print("PDB file read successfully. Starting search...")

#             async with aiohttp.ClientSession() as session:
#                 async with session.post('https://search.foldseek.com/api/ticket', data={
#                                 'q': data,
#                                 'database[]': ['afdb50', 'afdb-swissprot', 'pdb100', 'afdb-proteome', 'mgnify_esm30', 'gmgcl_id'],
#                                 'mode': '3diaa',
#                             }).json() as ticket_response:
#                     response_text = await ticket_response.text()
#                     print(f"Response from ticket request: {response_text}")
#                     if ticket_response.status != 200:
#                         raise ClientResponseError(ticket_response.request_info, ticket_response.history, status=ticket_response.status)
#                     try:
#                         ticket = await ticket_response.json()
#                     except aiohttp.ContentTypeError:
#                         print("Ticket response is not in JSON format")
#                         return response_text
#                     print("Search started. Ticket ID:", ticket['id'])

#                     # Collecting results
#                     repeat = True
#                     while repeat:
#                         status = await fetch(session, 'https://search.foldseek.com/api/ticket/' + ticket['id']).json()
#                         if isinstance(status, str):
#                             print("Status is not JSON format:", status)
#                             sys.exit(0)
#                         if status['status'] == "ERROR":
#                             # Handle error
#                             print("Error in status:", status)
#                             sys.exit(0)
#                         print("Status:", status)

#                         # Wait a short time between poll requests
#                         await asyncio.sleep(1)
#                         repeat = status['status'] != "COMPLETE"
#                     print("Search complete. Retrieving results...")

#                     # Get all hits for the first query (0)
#                     for _ in range(retries):
#                         try:
#                             result = await fetch(session, 'https://search.foldseek.com/api/result/' + ticket['id'] + '/0').json()
#                             if isinstance(result, str):
#                                 print("Result is not JSON format:", result)
#                                 continue
#                             print("Results retrieved successfully.")
#                             break
#                         except ClientResponseError as e:
#                             print(f"Error retrieving results: {e}. Retrying...")
#                             await asyncio.sleep(5)
#                     else:
#                         raise Exception("Failed to retrieve results after multiple attempts.")

#                     # Iterating over results
#                     foldseek_results_rows = []
#                     for results_dbs in result['results']:
#                         for alignment in results_dbs['alignments']:
#                             foldseek_results_rows.append(pd.DataFrame.from_dict({'pdb_file': [pdb_path],
#                                                                                 'db': [results_dbs['db']], 
#                                                                                 'query': [alignment.get('query', '-')],
#                                                                                 'target': [alignment.get('target', '-')],
#                                                                                 'seqId': [alignment.get('seqId', '-')],
#                                                                                 'alnLength': [alignment.get('alnLength', '-')],
#                                                                                 'missmatches': [alignment.get('missmatches', '-')],
#                                                                                 'gapsopened': [alignment.get('gapsopened', '-')],
#                                                                                 'qStartPos': [alignment.get('qStartPos', '-')],
#                                                                                 'qEndPos': [alignment.get('qEndPos', '-')],
#                                                                                 'dbStartPos': [alignment.get('dbStartPos', '-')],
#                                                                                 'dbEndPos': [alignment.get('dbEndPos', '-')],
#                                                                                 'prob': [alignment.get('prob', '-')],
#                                                                                 'eval': [alignment.get('eval', '-')],
#                                                                                 'score': [alignment.get('score', '-')],
#                                                                                 'qLen': [alignment.get('qLen', '-')],
#                                                                                 'dbLen': [alignment.get('dbLen', '-')],
#                                                                                 'qAln': [alignment.get('qAln', '-')],
#                                                                                 'dbAln': [alignment.get('dbAln', '-')],
#                                                                                 'tCa': [alignment.get('tCa', '-')],
#                                                                                 'tSeq': [alignment.get('tSeq', '-')],
#                                                                                 'taxId': [alignment.get('taxId', '-')],
#                                                                                 'taxName': [alignment.get('taxName', '-')]}))
#                     return pd.concat(foldseek_results_rows)
#         else:
#             print("Output file already exists.")
#     except Exception as e:
#         print(e)

# # 测试函数
# if __name__ == "__main__":
#     loop = asyncio.get_event_loop()
#     result_df = loop.run_until_complete(foldseek_search(pdb_path, output_file))
#     if result_df is not None:
#         print(result_df.head())
#     else:
#         print("No results returned.")
        